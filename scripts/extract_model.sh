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
	mkdir $2
	for f in `ls $1`; do
		java -cp ./ BuildContactMapFromPDB $1/$f $ANGS $BOOL
	done
	mv $1/*.cm $2
}
                                                        


##############
#extract_model $HOME/Downloads/pdb40d $HOME/Downloads/contact_maps_pdb40d
make_contact_map $HOME/Downloads/pdb40d $HOME/Downloads/contact_maps_pdb40d

                                         