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





def computeF(n):
    # print(n)
    t0 = get_time()
    time.sleep(1)
    t1 = get_time()
    time_elapsed(t0["tms"], t1["tms"])



if __name__ == "__main__":
    n = 10000
    while n <= 10000000:
        computeF(n)
        n *= 10


# ============================= 8-proc ===================================================================
#	n =    10000, digits =   35660, time =    0.145sec., memory 0.0M, max memory 0M
#	n =   100000, digits =   456574, time =   0.533sec., memory 0.0M, max memory 0M
#	n =  1000000, digits =  5565709, time =  22.055sec., memory 2.0M, max memory 2M
#	n = 10000000, digits = 65657060, time = 564.378sec., memory 31.0M, max memory 31M
# ========================================================================================================
