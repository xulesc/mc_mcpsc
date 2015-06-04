#!/bin/sh

extract_model()
{
	echo "extracting models"
	mkdir $2
	for f in `ls $1`; do
		./$SCRIPT_DIR/extractModel.pl $1/$f $2/$f
	done
}


##############
extract_model $HOME/Downloads/pdb40d $HOME/Downloads/contact_maps_pdb40d

                                         