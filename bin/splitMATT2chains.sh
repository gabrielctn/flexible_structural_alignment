#!/bin/bash

#
# splitMATT2chains.sh v. 1.0 [2018-02-20]
# d.a.suplatov@belozersky.msu.ru
#

#THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND
#PLEASE CHECK MANUALLY THE CONSISTENCY OF THE GENERATED OUTPUT
#PLEASE REPORT ANY PROBLEMS TO D.A.SUPLATOV@BELOZERSKY.MSU.RU

if [ ! -n "$3" ];then
    echo ""
    echo "Splits MATT pdb alignment into independent pdb files"
    echo "Usage: $0 <parMATT_output.pdb> <parMATT_output.txt> <output_dir>"
    echo ""
    exit
fi

PDB_INPUT=$1
TXT_INPUT=$2
OUTDIR=$3
PREFIX="$OUTDIR/"

echo
echo "THIS SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND"
echo "PLEASE CHECK MANUALLY THE CONSISTENCY OF THE GENERATED OUTPUT"
echo "PLEASE REPORT ANY PROBLEMS TO D.A.SUPLATOV@BELOZERSKY.MSU.RU"
echo ""
echo "Starting the script ... "
echo "Input: parMATT's information file is  $TXT_INPUT"
echo "Input: parMATT's 3D-alignment file is $PDB_INPUT"
echo "Input: separated PDBs will be printed to folder $OUTDIR"

if [ -d "$OUTDIR" ];then
    echo "Error: Output directory $OUTDIR already exists"
    exit
fi

mkdir $OUTDIR

awk 'BEGIN {i=1;a=1;c2p=1;pdbid=1;}

/^(ATOM|SEQRES|HETATM)/ {
    if (c2p!=0) {
        print "Info: Finished reading the parMATT information file", $TXT_INPUT
        print "Info: Reading the parMATT 3D-alignment file", $PDB_INPUT
        c2p=0;
    }
}

// {
    if (c2p!=0) {                        
           if ($0 == "") {                
                if (c2p==1) {
                    print "Info: Reading the parMATT information file", $TXT_INPUT
                }
                c2p=c2p+1; pdbid=1; next;
           }
           
           if (c2p==2) {
                PDBNAME=$1
                gsub(":", "_", PDBNAME)
                ID2PDB[pdbid]=PDBNAME
                print "Info: Chain", pdbid,"->", PDBNAME
                pdbid=pdbid+1
           }
    } 
}

/^ATOM/ {
    name="'$PREFIX'"ID2PDB[i]".pdb"; 
    na=sprintf("ATOM%7s", a)
    sub(/ATOM\s+/, "");
    sub($1, na)
    a=a+1
    print >> name    
}
/^TER/ {
    print "Info: Chain", i, "printed to file", name;
    i=i+1;
    a=1;
} ' $TXT_INPUT $PDB_INPUT

echo "Info: Finished reading the parMATT 3D-alignment file $PDB_INPUT"
echo "Info: The PDB files were printed to the folder $OUTDIR"
echo "Done!"
echo ""