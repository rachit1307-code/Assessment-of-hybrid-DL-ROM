#!/bin/bash

# nMOR Package directory and base directory
nMOR_DIR=/home/mechuser/Documents/MLforCFD/Github/OMAE_2022/Session_3/Workshop_3.2/nMOR
export PYTHONPATH=${nMOR_DIR}
BASE_DIR=`pwd`

# Model and data directories
OUTPUT_DIR=${BASE_DIR}/32P
DATA_DIR=${BASE_DIR}
NUM_UNITS=32
NUM_INFER_STEPS=1000
CKPT=500000

# Checkpoint directory
CKPT_DIR=/home/mechuser/Documents/MLforCFD/Github/OMAE_2022/Session_3/Workshop_3.2/32P/ckpt/nMOR_model.ckpt-${CKPT}

# Dataset to perform inference on and output directory
INFER_DATA_DIR=${DATA_DIR}/pres_rb_infer_2.h5
INFER_OUTPUT=pres_rb_pred_try.h5

python3 -m nMOR.nMOR \
    --num_units=${NUM_UNITS} \
    --num_infer_steps=${NUM_INFER_STEPS} \
    --time_major=True \
    --optimizer=adam \
    --out_dir=${OUTPUT_DIR} \
    --infer_data_file=${INFER_DATA_DIR} \
    --infer_output_file=${INFER_OUTPUT} \
    --data_size 64 64 \
    --batch_size=1 \
    --random_seed=1 \
    --num_gpus=0 \
    --ckpt=${CKPT_DIR} \
    --override_loaded_hparams=True \
    --num_intra_threads=8 \
--num_inter_threads=0
