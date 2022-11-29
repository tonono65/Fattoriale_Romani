import time 


def get_time():
    # tempo = {}
    tempo = dict()
    tempo["tms"] = time.time()
    tempo["tms_local"] = time.localtime(tempo["tms"])
    tempo["tms_str"] = time.strftime("%d-%m-%Y, %H:%M:%S", tempo["tms_local"])
    return  tempo
    

def time_elapsed(time_start, time_end):
    ore = int((time_end-time_start) // 3600)
    minuti = int(((time_end - time_start) - ore * 3600) // 60)
    secondi = (time_end - time_start) - ore * 3600 - minuti * 60 
    
    '''
    to_print = "Tempo di calcolo: {ore} ore   {min} minuti  {sec} secondi".format(
        ore=str(ore), min=str(minuti), sec=str(secondi)    
    ) 
    '''
    to_print = f"Tempo di calcolo: {ore} ore   {minuti} minuti  {secondi} secondi"
     

    print(to_print)
    # print("Tempo di calcolo: " + str(ore) + " ore " + str(minuti) + " minuti " + str(secondi) + " secondi")




