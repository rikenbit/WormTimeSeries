#!/bin/bash
# args n1_28sample
mkdir -p data/$@ 

# 既存の距離行列のmkdir&コピー
# n1_28sampleにはまだ24個体のサンプルしかない状態で実行
cp -r output/WTS4/n1_28sample output/WTS4/n1_24sample_add3
cp -r output/WTS4/n1_28sample output/WTS4/n1_24sample_add8
cp -r output/WTS4/n1_28sample output/WTS4/n1_24sample_add20
cp -r output/WTS4/n1_28sample output/WTS4/n1_24sample_add25

# n1_28sample配下に28個体の距離行列が作成された後に実行
cp output/WTS4/n1_28sample/all/SBD_abs/Distance/SampleNumber_3.RData output/WTS4/n1_24sample_add3/all/SBD_abs/Distance/SampleNumber_3.RData
cp output/WTS4/n1_28sample/stimAfter/SBD_abs/Distance/SampleNumber_3.RData output/WTS4/n1_24sample_add3/stimAfter/SBD_abs/Distance/SampleNumber_3.RData
cp output/WTS4/n1_28sample/all/EUCL/Distance/SampleNumber_3.RData output/WTS4/n1_24sample_add3/all/EUCL/Distance/SampleNumber_3.RData
cp output/WTS4/n1_28sample/stimAfter/EUCL/Distance/SampleNumber_3.RData output/WTS4/n1_24sample_add3/stimAfter/EUCL/Distance/SampleNumber_3.RData

cp output/WTS4/n1_28sample/all/SBD_abs/Distance/SampleNumber_8.RData output/WTS4/n1_24sample_add8/all/SBD_abs/Distance/SampleNumber_8.RData
cp output/WTS4/n1_28sample/stimAfter/SBD_abs/Distance/SampleNumber_8.RData output/WTS4/n1_24sample_add8/stimAfter/SBD_abs/Distance/SampleNumber_8.RData
cp output/WTS4/n1_28sample/all/EUCL/Distance/SampleNumber_8.RData output/WTS4/n1_24sample_add8/all/EUCL/Distance/SampleNumber_8.RData
cp output/WTS4/n1_28sample/stimAfter/EUCL/Distance/SampleNumber_8.RData output/WTS4/n1_24sample_add8/stimAfter/EUCL/Distance/SampleNumber_8.RData

cp output/WTS4/n1_28sample/all/SBD_abs/Distance/SampleNumber_20.RData output/WTS4/n1_24sample_add20/all/SBD_abs/Distance/SampleNumber_20.RData
cp output/WTS4/n1_28sample/stimAfter/SBD_abs/Distance/SampleNumber_20.RData output/WTS4/n1_24sample_add20/stimAfter/SBD_abs/Distance/SampleNumber_20.RData
cp output/WTS4/n1_28sample/all/EUCL/Distance/SampleNumber_20.RData output/WTS4/n1_24sample_add20/all/EUCL/Distance/SampleNumber_20.RData
cp output/WTS4/n1_28sample/stimAfter/EUCL/Distance/SampleNumber_20.RData output/WTS4/n1_24sample_add20/stimAfter/EUCL/Distance/SampleNumber_20.RData

cp output/WTS4/n1_28sample/all/SBD_abs/Distance/SampleNumber_25.RData output/WTS4/n1_24sample_add25/all/SBD_abs/Distance/SampleNumber_25.RData
cp output/WTS4/n1_28sample/stimAfter/SBD_abs/Distance/SampleNumber_25.RData output/WTS4/n1_24sample_add25/stimAfter/SBD_abs/Distance/SampleNumber_25.RData
cp output/WTS4/n1_28sample/all/EUCL/Distance/SampleNumber_25.RData output/WTS4/n1_24sample_add25/all/EUCL/Distance/SampleNumber_25.RData
cp output/WTS4/n1_28sample/stimAfter/EUCL/Distance/SampleNumber_25.RData output/WTS4/n1_24sample_add25/stimAfter/EUCL/Distance/SampleNumber_25.RData