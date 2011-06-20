#!/usr/bin/env python2.7
import threader as sc
from subprocess import call


if __name__ == '__main__':
  
  #for j in [10,15,20,25,30,35,40,45,50]:
  for k in range(1,15):
    j = 51
    q = ''
    for i in range(j):
      q += 'Q'
    #print len(q), q
    name = "PolyQ_" + str(j-1) + "_" + str(k)  + '.pdb'
    template = sc.MakeVariableLengthTemplate(j, 1)
  
    task2 = 'loopy -prm 2 -ini 10 -obj 1-' + str(j-1) + ' -cid A -out 14 -res ' + q + ' ' + template
    #print task2
    call(task2.split())
    nametemp = ((template.split('.pdb'))[0]) + '_loopy' + '.' + str(k) + '.pdb'
    task3 = "tail -n " + str((j-1)*17) + " " + nametemp + " > " + name + "temp"
    call(task3, shell=True)

    cuttask = "cut -c0-55 " + name + "temp" + " > " + name
    call(cuttask,shell=True)

    rmtask = "rm " + name + "temp " + nametemp
    call(rmtask.split())
    #rename = 'mv ' + nametemp + ' ' + name
    #print rename, '\n'
    #call(rename.split())


