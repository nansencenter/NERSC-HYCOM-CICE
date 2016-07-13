if [ $# -ne 1 ] ; then
   echo "This routine will back up one experiment to migrate"
   echo "You need to specify the experiment number "
   echo ""
   echo "Example:"
   echo "  $(basename $0) 01.2"
   echo "where 01.2 is the experiment you want to back up and"
   exit
fi





X=$1
source expt_${X}/EXPT.src
R=`pwd | awk --field-separator=/  '{print $(NF)}'`
E=$(echo $X | tr -d .)
user=`whoami`


echo ${E}
echo ${R}

#first delete file in SCRATCH folder that you do not want to save
rm -rf expt_${X}/SCRATCH/*

STRING="topo/partit/depth* topo/depth* topo/gr* topo/reg* Build_V${V}_X${X} REGION.src expt_${X}"

#here we ensure that the following exist before trying to tar them
[ -e topo/lat*                   ] &&  STRING=`echo $STRING topo/lat*` 
[ -e topo/newpos*                ] &&  STRING=`echo $STRING topo/newpos*` 
[ -e force/rivers/${E}/rivers.a  ] &&  STRING=`echo $STRING force/rivers/${E}/rivers.*` 
[ -e force/rivers/${E}/rivnitr.a ] &&  STRING=`echo $STRING force/rivers/${E}/rivnitr.*` 
[ -e force/rivers/${E}/rivphos.a ] &&  STRING=`echo $STRING force/rivers/${E}/rivphos.*` 
[ -e force/rivers/${E}/rivsili.a ] &&  STRING=`echo $STRING force/rivers/${E}/rivsili.*` 
[ -e force/seawifs/kpar.a        ] &&  STRING=`echo $STRING force/seawifs/kpar.*`
[ -e topo/tbaric.a               ] &&  STRING=`echo $STRING topo/tbaric.*`
[ -e  tides_nersc/${E}           ] &&  STRING=`echo $STRING tides_nersc/${E}/*`
[ -e  relax/${E}                 ] &&  STRING=`echo $STRING relax/${E}/relax_*`
[ -e  nest_nersc/${E}/outer      ] &&  STRING=`echo $STRING nest_nersc/${E}/outer/*`
[ -e  nest_nersc/${E}/inner      ] &&  STRING=`echo $STRING nest_nersc/${E}/inner/*`



#tar the selected file
tar -cvf backup_expt${X}.tar    $STRING || { echo "problem during the tar we stop";exit 1;  }
gzip  backup_expt${X}.tar || { echo "problem during the zipping we stop";exit 1;  }
if ! [ -w ${BACKUP_PATH} ]; then echo "The BACKUP_PATH set in the
   REGION.src does not exist, I quit";exit 1;fi
mkdir -p ${BACKUP_PATH}/${R}/
mv  backup_expt${X}.tar.gz  ${BACKUP_PATH}/${R}/ || { echo "problem during the moving we stop";exit 1; }
rm -rf backup_expt${X}.tar.gz || { echo "problem during the deleting";exit 1; }
