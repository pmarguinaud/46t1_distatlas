#!/bin/bash

set -x

~/SAVE/mpiauto/mpiauto --wrap --wrap-stdeo -np 4 -- ./bin/ATLAS_ARPEGE_F 
