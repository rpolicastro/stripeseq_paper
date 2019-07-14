#!/bin/bash

CONTAINER="library://rpolicastro/default/stripe_seq:1.0.0"

##########################################
## Use Singularity Container for Analysis
##########################################

## Check if container is downloaded.

CONTAINER_NAME=$(basename $CONTAINER | tr ":" "_").sif

if [ ! -f "singularity/${CONTAINER_NAME}" ]; then
	singularity pull $CONTAINER 
	mv $CONTAINER_NAME ./singularity/$CONTAINER_NAME
fi

## Start container.

singularity shell \
-eCB $PWD -H $PWD \
singularity/$CONTAINER_NAME
