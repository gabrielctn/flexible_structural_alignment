#!/usr/bin/perl -w

# Copyright Jean-Christophe Gelly (April 10 2013)
#
# jean-christophe.gelly@univ-paris-diderot.fr
#
# This software is a computer program whose purpose is to 
# compute several measures related to structural alignment quality
# from a given set of aligned protein structure in PDB format.
#
# This software is governed by the CeCILL-B license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL-B
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL-B license and that you accept its terms.

# GDT v0.1 2012 04 10

use strict;
### 1) PARSE ARGVS ####
#---------------------#
my $argv=join(' ',@ARGV);
chomp $argv;

my @tab_tmp=split(' ',$argv);

my @tab_pdb;
foreach my $tmp (@tab_tmp)
{
    if ($tmp =~/ / or $tmp =~/^$/)
    {
    }
    else
    {
        push @tab_pdb,$tmp;
    }
}

##### 2) Extract DATA FROM PDB FILES ######
#-----------------------------------------#
my @tab_seq;
my @tab_CA;
my @tab_model;
my $num=0;
for (my $i=0 ; $i <= $#tab_pdb; $i++)
{
    my $model=0;
    my $pdb="$tab_pdb[$i]";
    my $name_pdb="$tab_pdb[$i]"."_model_"."$model";
    my $seq="";
    my @tab_coo;

    #PARSE ATOM
    my $check_ter=0;
    open(F,"$pdb") or die "Cannot open $pdb: $!\n";
    #Check if there are more than one model
    WHILE:while(my $line=<F>)
	{
		my $record_name = substr($line,  0, 6);  # columns 1-6 = ATOM or TER
        #Check if there are more than one model
        if ($record_name =~ /^TER/)
        {
            $tab_CA[$num]=[@tab_coo];
            $tab_seq[$num]=$seq;
            $tab_model[$num]=$name_pdb;
            $num++;
            $model++;
            $name_pdb="$pdb"."_model_"."$model";
            $seq="";
            @tab_coo=();
            $check_ter=1;
        }
        #Analyse only ATOM fields
	    if ($record_name =~/^ATOM/)
        {
		    my $number      = substr($line,  6, 5);  # columns 7-11
		    my $atom_name   = substr($line, 12, 3);  # columns 13-16
		    my $alt_loc     = substr($line, 16, 1);  # columns 17

            #Erase space
		    $atom_name=~s/\s//g;

            #Skip alternate position
		    if ($alt_loc !~ / /){
			    next FOR if ($alt_loc !~/^1|A$/);
		    }
		    my $aa          = substr($line, 17, 3);  # columns 18-20
		    my $chain       = substr($line, 21, 1);  # columns 22
		    my $no_aa       = substr($line, 22, 4);  # columns 23-26
		    my $insert_code = substr($line, 26, 4);  # columns 27
		    my $x           = substr($line, 30, 8);  # columns 31-38
		    my $y           = substr($line, 38, 8);  # columns 39-46
		    my $z           = substr($line, 46, 8);  # columns 47-54

            $x=~s/\s+//g;
            $y=~s/\s+//g;
            $z=~s/\s+//g;


            next WHILE if ($atom_name !~/^CA$/);


            #Keep coordinates
            push @tab_coo,[("$x", "$y" ,"$z")];
            # INSERT SEQ
            my $aa1=three2one($aa);
            $seq.=$aa1;
	    }
    }
    if ($check_ter == 0)
    {
        $tab_CA[$num]=[@tab_coo];
        $tab_seq[$num]=$seq;
        $tab_model[$num]=$name_pdb;
        $num++;
    }

}

###### MAIN ######
#----------------#

my $number=1;
for (my $i=0 ; $i < $#tab_model; $i++)
{
	my $pdb1=$tab_model[$i];
#   print "Chain 1 : $pdb1\n";
    my $seq1=$tab_seq[$i];
	for (my $j=$i+1 ; $j <= $#tab_model; $j++)
	{
        printf ("#### Alignment %3d ####\n",$number);
        my $pdb2=$tab_model[$j];
        my $seq2=$tab_seq[$j];
        my $ref_tab_tab_distance           = all_distance($tab_CA[$i],$tab_CA[$j],$tab_seq[$i],$tab_seq[$j]);
        my ($ref_tab_dist,$align1,$align2) = dp($seq1,$seq2,$pdb1,$pdb2,$ref_tab_tab_distance);
        $number++;
    }
}
#print STDERR "ok\n";


###### COMPUTE ALL DISTANCE BETWEEN ATOMS  #########
#--------------------------------------------------#
sub all_distance
{
    my $ref_tab_atom1=shift;
    my $ref_tab_atom2=shift;


    my $seq1=shift;
    my $seq2=shift;

    my @tab_seq1=split('',$seq1);
    my @tab_seq2=split('',$seq2);


    my @tab_tab_distance;

    for (my $i=0 ; $i <= $#$ref_tab_atom1; $i++)
    {
        my $ref_coo_atom1=$$ref_tab_atom1[$i];
        my $x1=$$ref_coo_atom1[0];
        my $y1=$$ref_coo_atom1[1];
        my $z1=$$ref_coo_atom1[2];
        my @tab_distance1;

#    printf("%-5s ",$tab_seq1[$i]);
        for (my $j=0 ; $j <= $#$ref_tab_atom2 ; $j++)
        {
            my $ref_coo_atom2=$$ref_tab_atom2[$j];
            my $x2=$$ref_coo_atom2[0];
            my $y2=$$ref_coo_atom2[1];
            my $z2=$$ref_coo_atom2[2];
            #print "$x1,$y1,$z1,$x2,$y2,$z2\n";
            my $distance=distance($x1,$y1,$z1,$x2,$y2,$z2);
#           printf ("%5.1f ",$distance);
            push @tab_distance1,$distance;
        }
        $tab_tab_distance[$i] = [@tab_distance1];
#        print "\n";
    }
    return(\@tab_tab_distance)
}


###### COMPUTE DISTANCE BETWEEN ATOM  #########
#---------------------------------------------#
sub distance
{
	my $x1=shift;
	my $y1=shift;
	my $z1=shift;
	my $x2=shift;
	my $y2=shift;
	my $z2=shift;
	my $dx2 = ($x1-$x2)**2;
	my $dy2 = ($y1-$y2)**2;
	my $dz2 = ($z1-$z2)**2;
	my $dt2 = $dx2+$dy2+$dz2;
	my $dt  = sqrt ($dt2);
	return($dt);
}


###### DYNAMIC PROGRAMMING AND COMPUTE MEASURES #########
#-------------------------------------------------------#
sub dp
{

	my $seq1=shift;
	my $seq2=shift;
    my $pdb1=shift;
    my $pdb2=shift;
	my $ref_tab_tab_distance=shift;
	my @tab_threshold   =qw(1 2 4 8);
    my @tab_threshold_HA=qw(0.5 1 2 4);

# Scoring scheme

	my $match    =  0; # Inititialized to given value

# Initialization

	my @matrix;

    $matrix[0][0]{score}{tmscore}   = 0;

	$matrix[0][0]{pointer} = 0; # -1

	for(my $j = 1; $j <= (length($seq2)+1); $j++) 
    {
        $matrix[0][$j]{score}{tmscore} = 0;
		$matrix[0][$j]{pointer} = 1; # 1 "left"
	}
	for (my $i = 1; $i <= (length($seq1)+1); $i++) 
    {
        $matrix[$i][0]{score}{tmscore}    = 0;
		$matrix[$i][0]{pointer}      = 3; # 3 "up"
	}

# fill

    my $size1=length($seq2);
    my $size2=length($seq1);

    #DEFINE D0 (from TM-align paper and TM-align source code (Version 20120707)
    my $size=$size2;

    #Normalize by smallest seq
    if ($size1 < $size2)
    {
        $size=$size1;
    }

    my $d0=(1.24*($size-15)**(1/3))-1.8;

    my $true_d0=$d0;
    if ($d0 < 0.5) 
    {
        $d0=0.5;
    }

    $d0=$d0+0.8; #Best for search from TM-align source code (Version 20120707)

	for (my $i = 1; $i <= length($seq1); $i++)
	{
		for (my $j = 1; $j <= length($seq2); $j++)
		{
			my ($diagonal_score, $left_score, $up_score);

            # Calculate match score
            #----------------------

            $match=${${$ref_tab_tab_distance}[$i-1]}[$j-1];

			$diagonal_score = 0;

            $diagonal_score = 1 / $size * ( $matrix[$i-1][$j-1]{score}{tmscore} + ( 1/(1+(($match/$d0)**2)) ));

# Calculate gap scores
#---------------------

#### UP SCORE
			$up_score=0;

            $up_score=1/$size*( $matrix[$i-1][$j]{score}{tmscore});

#### LEFT SCORE
			$left_score=0;

            $left_score=1/$size*( $matrix[$i][$j-1]{score}{tmscore});

# Choose best score
#-------------------

			if ($diagonal_score >= $up_score and $diagonal_score >= $left_score) 
			{

                $matrix[$i][$j]{score}{tmscore}= $matrix[$i-1][$j-1]{score}{tmscore} + ( 1/(1+($match/$d0)**2) );

				$matrix[$i][$j]{pointer} = 2; # 2 "diagonal"
			}
			elsif ($up_score <= $left_score)
			{
				$matrix[$i][$j]{pointer} = 1; # 1 left

                $matrix[$i][$j]{score}{tmscore}=$matrix[$i][$j-1]{score}{tmscore};

			}
			else
			{

                $matrix[$i][$j]{score}{tmscore}=$matrix[$i-1][$j]{score}{tmscore};

				$matrix[$i][$j]{pointer} = 3; # 3 up
			}
		}
	}


# Trace-back

	my $align1 = "";
	my $align2 = "";

# Start at last cell of matrix
	my $i = length($seq1);
	my $j = length($seq2);

	my $align3="";

    printf ("Chain 1: $pdb1\n");
    printf ("Chain 2: $pdb2\n");

    my $tmscore_search=1/$size *( $matrix[$i][$j]{score}{tmscore});
    printf ("TM-score-Search (Normalized on Smallest chain): %5.3f\n",$tmscore_search);


    my @tab_threshold_HR;
    for (my $i=0 ; $i <= $#tab_threshold ; $i++)
    {
        my $step=($tab_threshold[$i]-$tab_threshold[$i-1])/100;
        for (my $j=$tab_threshold[$i-1]; $j <= $tab_threshold[$i]; $j=$j+$step)
        {
            push @tab_threshold_HR,$j;
        }
    }

    my $gdt_HR=0;
    my $gdt_HA=0;
    my $gdt=0;
    my $rmsd=0;
    my $num_ali=0;
    my $num_ali_rmsd=0;
    my $tmscore1=0;
    my $tmscore2=0;
    my $pid=0;
    my $pid_rmsd=0;

    #Compute d0
    my $d0_1=(1.24*($size1-15)**(1/3))-1.8;
    my $d0_2=(1.24*($size2-15)**(1/3))-1.8;

    if ($d0_1 < 0.5) 
    {
        $d0_1=0.5;
    }

    if ($d0_2 < 0.5) 
    {
        $d0_2=0.5;
    }


	while (1) 
	{
		last if ($matrix[$i][$j]{pointer} == 0); # ends at first cell of matrix
		if ($matrix[$i][$j]{pointer} == 2)
		{
			my $aa1=substr($seq1, $i-1, 1);
			my $aa2=substr($seq2, $j-1, 1);
			my $dist=${${$ref_tab_tab_distance}[$i-1]}[$j-1];
			my $score=sprintf("%1d",$dist);
			if ($dist > 5.0) {$aa1=lc($aa1);$aa2=lc($aa2);}
            if ($dist >=10.0) {$score="x";}
			$align1 .= $aa1;
			$align2 .= $aa2;
			$align3 .= $score;

            #Number aligned
            $num_ali++;

            #PID
            if ($aa1 eq $aa2)
            {
                $pid++;
            }

            #RMSD
            my $dist_rmsd=$dist;
            if($dist <= 5)
            {
                $rmsd=$rmsd+(($dist_rmsd)**2);
                $num_ali_rmsd++;
                if ($aa1 eq $aa2){$pid_rmsd++;}
            }
            #Classical GDT
	        for (my $a=0 ; $a <= $#tab_threshold ; $a++) 
            {
                if ($dist <= $tab_threshold[$a])
                {
                    $gdt++;
                }
            }
            #HA GDT
	        for (my $a=0 ; $a <= $#tab_threshold_HA ; $a++) 
            {
                if ($dist <= $tab_threshold_HA[$a])
                {
                    $gdt_HA++;
                }
            }

            #HR GDT
	        for (my $a=0 ; $a <= $#tab_threshold_HR ; $a++) 
            {
                if ($dist <= $tab_threshold_HR[$a])
                {
                    $gdt_HR++;
                }
            }

            #TM-SCORE
            $tmscore1=$tmscore1 + ( 1/(1+(($dist/$d0_1)**2)));
            $tmscore2=$tmscore2 + ( 1/(1+(($dist/$d0_2)**2)));

			$i--;
			$j--;
		}
		elsif ($matrix[$i][$j]{pointer} == 3)
		{
			$align1 .= substr($seq1, $i-1, 1);
			$align2 .= "-";
			$align3 .= " ";
			$i--;
		}
		elsif ($matrix[$i][$j]{pointer} == 1)
		{
			$align1 .= "-";
			$align2 .= substr($seq2, $j-1, 1);
			$align3 .= " ";
			$j--;
		}
	}

    #Number of aa
    printf("Number of AAs Sequence 1  : %d\n",length($seq1));
    printf("Number of AAs Sequence 2  : %d\n",length($seq2));


    #Number aligned
    printf("Number of residues aligned: %d\n",$num_ali);
    printf("Number of residues aligned: %d (below a threshod of 5.0 A)\n",$num_ali_rmsd);

    # Percentage of identity
    printf("Percentage of identity    : %5.2f\n",$pid/$num_ali*100);
    printf("Percentage of identity    : %5.2f (on residues below a threshod of 5.0 A)\n",$pid_rmsd/$num_ali_rmsd*100);

    #Compute RMSD
    $rmsd=sqrt($rmsd/$num_ali_rmsd);
    printf("RMSD     : %5.3f (on residues below a threshod of 5.0 A)\n",$rmsd);

    #Compute TM-SCORE
    $tmscore1=1/$size1 *  $tmscore1;
    $tmscore2=1/$size2 *  $tmscore2;
    printf("TM-score : %5.3f (if normalized by length of Chain 1)\n",$tmscore1);
    printf("TM-score : %5.3f (if normalized by length of Chain 2)\n",$tmscore2);


    #Compute GDT (Normalized by smallest chain
    $gdt   = 100 * ($gdt/(($#tab_threshold+1)*$size));
    $gdt_HA= 100 * ($gdt_HA/(($#tab_threshold_HA+1)*$size));
    $gdt_HR= 100 * ($gdt_HR/(($#tab_threshold_HR+1)*$size));
    printf("GDT      : %5.1f (Normalized by smallest chain)\n",$gdt);
    printf("GDT-HA   : %5.1f (Normalized by smallest chain)\n",$gdt_HA);
    printf("GDT-HR   : %5.1f (Normalized by smallest chain)\n",$gdt_HR);


	$align1 = reverse $align1;
	$align2 = reverse $align2;
	$align3 = reverse $align3;
	print "CHAIN_1:$align1\n";
	print "CHAIN_2:$align2\n";
	print "DIST   :$align3\n";
    print "\n";
	return($align1,$align2);
}

###### CONVERT AA 3to1 #########
#------------------------------#
sub three2one 
{ 
# converts a 3-letter-code to a 1-letter-abbreviation

	my $Code=shift;
	$Code =~ s/\s//g;  # remove all blanks
#Verify if lowercase
	my $check_lowercase=0;
	if ($Code =~/[a-z]{3}/)
	{
		$check_lowercase=1;
	}
	$Code = uc ($Code);
#    sleep 1;
#    print "CODE|$Code|\n";
#    sleep 1;

	if    ($Code eq 'ALA') {$Code = 'A'} # Alanine
	elsif ($Code eq 'ASX') {$Code = 'B'} # Aspartic acid or Asparagine
	elsif ($Code eq 'CYS') {$Code = 'C'} # Cysteine
	elsif ($Code eq 'ASP') {$Code = 'D'} # Aspartic acid
	elsif ($Code eq 'GLU') {$Code = 'E'} # Glutamic acid
	elsif ($Code eq 'PHE') {$Code = 'F'} # Phenylalanine
	elsif ($Code eq 'GLY') {$Code = 'G'} # Glycine
	elsif ($Code eq 'HIS') {$Code = 'H'} # Histidine
	elsif ($Code eq 'ILE') {$Code = 'I'} # Isoleucine
	elsif ($Code eq 'LYS') {$Code = 'K'} # Lysine
	elsif ($Code eq 'LEU') {$Code = 'L'} # Leucine
	elsif ($Code eq 'MET') {$Code = 'M'} # Methionine
	elsif ($Code eq 'ASN') {$Code = 'N'} # Asparagine
	elsif ($Code eq 'PRO') {$Code = 'P'} # Proline
	elsif ($Code eq 'GLN') {$Code = 'Q'} # Glutamine
	elsif ($Code eq 'ARG') {$Code = 'R'} # Arginine
	elsif ($Code eq 'SER') {$Code = 'S'} # Serine
	elsif ($Code eq 'THR') {$Code = 'T'} # Threonine
	elsif ($Code eq 'SEC') {$Code = 'U'} # Selenocysteine
	elsif ($Code eq 'VAL') {$Code = 'V'} # Valine
	elsif ($Code eq 'TRP') {$Code = 'W'} # Tryptophan
	elsif ($Code eq 'UNK') {$Code = 'X'} # Unknown amino acid
	elsif ($Code eq 'TYR') {$Code = 'Y'} # Tyrosine
	elsif ($Code eq 'GLX') {$Code = 'Z'} # Glutamic acid or Glutamine
	else 
	{
		return ("X");
	}
	if ($check_lowercase==1) {$Code=lc($Code);};
	return($Code);
}


