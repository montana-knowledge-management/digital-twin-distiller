import pandas as pd
from digital_twin_distiller import ModelDir

cogging = pd.read_csv( ModelDir.DATA / "df_cogging.csv")
locked