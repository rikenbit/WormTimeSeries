source("src/functions_WTS5_ts2net_dist.R")

#### args setting####
#### test args####
args_sample <- c("2")
args_neuron <- c("data/normalize_1/ReadData_2.RData")
args_stim_xlsx <- c("data/stimulation/stimulation_timing.xlsx")

load(args_neuron)
eval(parse(text=paste0("ReadData <- ReadData_",args_sample)))
ReadData <- .ReadData_stimAfter(ReadData, args_stim_xlsx)

#### tssim_event_sync関数####
# https://cran.r-project.org/web/packages/ts2net/ts2net.pdf
# minimaldata

# test_tts1 <- ReadData[,1]
t_tts1 <- ReadData["ASER"]$ASER
t_tts2 <- ReadData["ASEL"]$ASEL

system.time(dist_value <- ts2net:::tssim_event_sync(t_tts1, t_tts2, tau_max = 500))

# 500TFで1秒@M1 Mac
t_tts1_500 <- t_tts1[1:500]
t_tts2_500 <- t_tts2[1:500]
system.time(dist_value <- ts2net:::tssim_event_sync(t_tts1_500, t_tts2_500, tau_max = 500))
  # ユーザ   システム       経過
  #    0.992      0.004      0.996

# 2000TFで16秒@M1 Mac
t_tts1_2000 <- t_tts1[1:2000]
t_tts2_2000 <- t_tts2[1:2000]
system.time(dist_value <- ts2net:::tssim_event_sync(t_tts1_2000, t_tts2_2000, tau_max = 500))
# ユーザ   システム       経過  
# 16.649      0.043     16.690 

# 6000TFで3分@M1 Mac
t_tts1 <- ReadData["ASER"]$ASER
t_tts2 <- ReadData["ASEL"]$ASEL

system.time(dist_value <- ts2net:::tssim_event_sync(t_tts1, t_tts2, tau_max = 500))
# ユーザ   システム       経過  
# 183.119      0.357    183.404 

# 200cell x 200cell 40000 *3分 =120000分=2000時間



#### tsdist_nmi関数 テスト####
# https://cran.r-project.org/web/packages/ts2net/ts2net.pdf
# test_tts1 <- ReadData[,1]
t_tts1 <- ReadData["ASER"]$ASER
t_tts2 <- ReadData["ASEL"]$ASEL

system.time(dist_value <- ts2net:::tsdist_nmi(t_tts1, t_tts2))
