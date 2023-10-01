# *********************************************************
#                                                     
# shell script to do the following:
# 
# - replace $(CC) with mpic++ in pluto the makefile
# - add $(CLOUDY_OBJ) to the objects being linked to pluto
#
# 2013 May 24th - M.Salz
#
# *********************************************************

sed -i ':a;N;$!ba;s/pluto: $(OBJ) \n\t$(CC) $(OBJ)/pluto: $(OBJ) \n\tmpic++ $(OBJ) $(CLOUDY_OBJ)/g' makefile

