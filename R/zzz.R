# The DIC helpers in cal_DIC_v7.R read the data and its dimensions as free
# variables; MINDS_algorithm() supplies them by rebinding the helpers into a
# local `dic_env` holding this call's data. Declare them here so R CMD check
# does not flag them as undefined globals.
globalVariables(c(
  "M_1", "M_2", "Nb", "Nc", "Nd1", "Nd2", "y_1", "y_2"
))
