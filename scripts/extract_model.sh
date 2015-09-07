#!/bin/sh

## Contact Map Parameters
ANGS='6.5'                      ## as specified in the paper
BOOL='false'
##

extract_model()
{
	echo "extracting models"
	mkdir $2
	for f in `ls $1`; do
		./$SCRIPT_DIR/extractModel.pl $1/$f $2/$f
	done
}

make_contact_map()
{
	echo "making contact maps"
	[ -e $2 ] || mkdir $2
	for f in `ls $1`; do
		java -cp ./ BuildContactMapFromPDB $1/$f $ANGS $BOOL
	done
	mv $1/*.cm $2
}
                                                        


##############
# extract_model $HOME/Downloads/pdb40d_full_dom $HOME/Downloads/pdb40d
#make_contact_map $HOME/Downloads/pdb40d $HOME/Downloads/contact_maps_pdb40d
#extract_model $HOME/Downloads/astral_40p $HOME/Downloads/astral_40p_models
#make_contact_map $HOME/Downloads/astral_40p_models $HOME/Downloads/contact_maps_astral_40p
#extract_model /media/lvm/astral/pdb_style/pdb40d /media/lvm/astral/pdb_style/pdb40d_models
#make_contact_map /media/lvm/astral/pdb_style/pdb40d_models /media/lvm/astral/pdb_style/pdb40d/pdb40d_models_contact_maps
#extract_model /media/lvm/astral/pdb_style/apdb40 /media/lvm/astral/pdb_style/apdb40_models
make_contact_map /media/lvm/astral/pdb_style/apdb40_models /media/lvm/astral/pdb_style/apdb40_models_contact_maps


                                         