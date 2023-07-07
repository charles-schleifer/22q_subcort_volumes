#!/bin/bash

# setup freesurfer
export FREESURFER_HOME=/Applications/freesurfer/
export SUBJECTS_DIR=$FREESURFER_HOME/subjects
source $FREESURFER_HOME/SetUpFreeSurfer.sh 

# mount hoffman cluster
#hoffman="/Users/charlie/Desktop/hoffman_mount"
#umount -f $hoffman
#echo "...mounting ${hoffman}" 
#sshfs schleife@hoffman2.idre.ucla.edu:/u/project/cbearden/data ${hoffman}

# study dir
#studydir="${hoffman}/22q_T1w_all/sessions_recon-all/"
studydir="/Users/charlie/Documents/qc_fs_subcort/"
seshpaths=$(ls -d ${studydir}/Q_*long*)

# slices to capture
xvals=(-25 -15 15 25)
yvals=(-25 -20 -15 -10 0)
zvals=(-20 0 5 10)

# set to 1 if want to copy from hoffman to local
copy=0
freeview=1

# loop through subjects in each study
for seshpath in $seshpaths; do
	sesh=$(basename $seshpath)
	echo $sesh
	# set up paths
	#qdcir=${studydir}/qc_fs_subcort/${sesh}
	qcdir=~/Documents/qc_fs_subcort/${sesh}
	qcfile=${qcdir}/freeview_commands.txt
	t1dir=${studydir}/${sesh}/mri
	t1mgz=T1.mgz
	thalmgz=ThalamicNuclei.long.mgz
	rhhipmgz=rh.hippoAmygLabels.long.mgz
	lhhipmgz=lh.hippoAmygLabels.long.mgz
	# copy files to local
	if [[  $copy -eq 1  ]]; then
		#cp -v ${t1dir}/${t1mgz} $qcdir
		#cp -v ${t1dir}/${thalmgz} $qcdir
		#cp -v ${t1dir}/${rhhipmgz} $qcdir
		#cp -v ${t1dir}/${lhhipmgz} $qcdir
		
		rsync -az schleife@hoffman2.idre.ucla.edu:/u/project/cbearden/data/22q_T1w_all/sessions_recon-all/${sesh}/mri/${t1mgz} $qcdir
		rsync -az schleife@hoffman2.idre.ucla.edu:/u/project/cbearden/data/22q_T1w_all/sessions_recon-all/${sesh}/mri/${thalmgz} $qcdir
		rsync -az schleife@hoffman2.idre.ucla.edu:/u/project/cbearden/data/22q_T1w_all/sessions_recon-all/${sesh}/mri/${rhhipmgz} $qcdir
		rsync -az schleife@hoffman2.idre.ucla.edu:/u/project/cbearden/data/22q_T1w_all/sessions_recon-all/${sesh}/mri/${lhhipmgz} $qcdir

	fi
	# create freeview command file
	if [[  $freeview -eq 1  ]]; then
		mkdir -p $qcdir
		echo "freeview -v ${qcdir}/mri/${t1mgz}" > $qcfile
		echo "freeview -v ${qcdir}/mri/${thalmgz}:colormap=LUT" >> $qcfile
		echo "freeview -v ${qcdir}/mri/${rhhipmgz}:colormap=LUT" >> $qcfile
		echo "freeview -v ${qcdir}/mri/${lhhipmgz}:colormap=LUT" >> $qcfile
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
	fi
done

