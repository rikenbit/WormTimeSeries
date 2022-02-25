#!/bin/bash
# args n1_28sample
mkdir -p data/$@ 

# 既存の距離行列のmkdir&コピー
# n1_28sampleにはまだ24個体のサンプルしかない状態で実行
cp -r output/WTS4/n1_28sample output/WTS4/n1_24sample_add3
cp -r output/WTS4/n1_28sample output/WTS4/n1_24sample_add8
cp -r output/WTS4/n1_28sample output/WTS4/n1_24sample_add20
cp -r output/WTS4/n1_28sample output/WTS4/n1_24sample_add25