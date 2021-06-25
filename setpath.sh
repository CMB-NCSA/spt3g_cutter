#!/usr/bin/env bash

PRODUCT_DIR="$1"

if [ $PRODUCT_DIR = "." ]; then
    PRODUCT_DIR=$PWD
fi

echo "Adding: $PRODUCT_DIR to paths"
export SPT3G_CUTTER_DIR=$PRODUCT_DIR
export PYTHONPATH=${PRODUCT_DIR}/python:${PYTHONPATH}
export PATH=${PRODUCT_DIR}/bin:${PATH}
