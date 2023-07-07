#!/bin/bash

# setup freesurfer
export FREESURFER_HOME=/Applications/freesurfer/
export SUBJECTS_DIR=$FREESURFER_HOME/subjects
source $FREESURFER_HOME/SetUpFreeSurfer.sh 

# mount hoffman cluster
hoffman="/Users/charlie/Desktop/hoffman_mount"
umount -f $hoffman
echo "...mounting ${hoffman}" 
sshfs schleife@hoffman2.idre.ucla.edu:/u/project/cbearden/data ${hoffman}

# study dirs
trio="/22q/qunex_studyfolder/sessions/"
prisma="/22qPrisma/qunex_studyfolder/sessions/"
triosmri="/22q/qunex_studyfolder/sessions_sMRIonly/"
prismasmri="/22qPrisma/qunex_studyfolder/sessions_sMRIonly/"


# slices to capture
xvals=(-30 -25 -15 -5 0 5 15 25 30)
yvals=(-30 -25 -20 -15 -10 0)
zvals=(-25 -20 -10 0 5 10)

# loop through subjects in each study
#for studydir in $trio $prisma $triosmri $prismasmri; do
for study in $triosmri $prismasmri; do
	case $study in
		${triosmri})
		studydir=${hoffman}/${triosmri}
		sessions=$(ls -d ${studydir}/Q_*)
		;;
		${prismasmri})
		studydir=${hoffman}/${prismasmri}
		sessions=$(ls -d ${studydir}/Q_*)
		;;		
		${trio})
		studydir=${hoffman}/${trio}
		sessions=$(ls -d ${studydir}/Q_*)
		;;
		${prisma})
		studydir=${hoffman}/${prisma}
		sessions=$(ls -d ${studydir}/Q_*)
		;;
	esac
	for seshpath in $sessions; do
		sesh=$(basename $seshpath)
		echo $sesh
		# set up paths
		#qdcir=${studydir}/qc_fs_subcort/${sesh}
		qcdir=~/Documents/qc_fs_subcort/${sesh}
		qcfile=${qcdir}/freeview_commands.txt
		t1dir=${studydir}/${sesh}/hcp/${sesh}/T1w/${sesh}/mri
		t1mgz=${t1dir}/T1.mgz
		thalmgz=${t1dir}/ThalamicNuclei.mgz
		rhhipmgz=${t1dir}/rh.hippoAmygLabels.mgz
		lhhipmgz=${t1dir}/lh.hippoAmygLabels.mgz
		# create freeview command file
		mkdir -p $qcdir
		echo "freeview -v ${t1mgz}" > $qcfile
		echo "freeview -v ${thalmgz}:colormap=LUT" >> $qcfile
		echo "freeview -v ${rhhipmgz}:colormap=LUT" >> $qcfile
		echo "freeview -v ${lhhipmgz}:colormap=LUT" >> $qcfile
		echo "-zoom 2.5" >> $qcfile
		echo "-viewport axial" >> $qcfile
		i=1
		for z in ${zvals[@]}; do 
			echo "-ras 0 0 ${z}" >> $qcfile
			echo "-ss ${qcdir}/axial${i}_z_${z}.png" >> $qcfile
			i=$(($i+1))
		done
		echo "-viewport sagittal" >> $qcfile
		i=1
		for x in ${xvals[@]}; do 
			echo "-ras ${x} 0 0" >> $qcfile
			echo "-ss ${qcdir}/sagittal${i}_x_${x}.png" >> $qcfile
			i=$(($i+1))
		done	
		echo "-viewport coronal" >> $qcfile
		i=1
		for y in ${yvals[@]}; do 
			echo "-ras 0 ${y} 0" >> $qcfile
			echo "-ss ${qcdir}/coronal${i}_y_${y}.png" >> $qcfile
			i=$(($i+1))
		done
		# get snapshots
		echo "-quit" >> $qcfile
		echo "...generating snapshots for ${sesh}"
		freeview -cmd $qcfile
	done
done
