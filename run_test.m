%rm=memlsmain(21, 50, 0, 0, 'test1layer.txt', 0, 273, 11)

%ra=amemlsmain(21, 50, 0, 0, 0, 0, 'test1layer.txt', 0, 273, 12, 0.1, 0.1)

rl=lmain2(10, 100, 50, 0, 0, 'test1layer.txt', 273, 11)

%rf=fmain2(35, 0, 80, 0, 0, 'test1layer.txt', 273, 11)

close all %close connection to graphic windows (or octave will crash)


run_MEMLS_Active_v2(10,20, 0.05,0.1,0)