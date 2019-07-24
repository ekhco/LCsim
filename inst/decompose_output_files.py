from __future__ import division
from ipdb import set_trace
import numpy as np
import os
from shutil import move
import socket
from sys import platform
import sys

def parse_results(year):
    global sims
    people = []
    people_ages = []
    people_cpd = []
    # 'output_' + str(year) + '.out'
    sims = open(sys.argv[1], 'r').read()
    if sims.find('<RUN>') != -1:
        sims = sims[sims.find('<RUN>') + 5:sims.find('</RUN>')].split('\n')
    else:
        sims = sims.split('\n')
    sims = [sim.split(';')[0:-1] for sim in sims if len(sim) != 0]
    for sim in sims:
        person = {}
        person['race'] = sim[0]
        person['sex'] = sim[1]
        person['yob'] = sim[2]
        person['init_age'] = sim[3]
        person['cess_age'] = sim[4]
        person['ocd_age'] = sim[5]
        if person['init_age'] != '-999':
            smoking_history = np.array(sim[6:])
            years = smoking_history.shape[0] // 2
            person['smoking_history'] = smoking_history.reshape(years, 2)
            people_ages.append(np.array(person['smoking_history'], dtype=np.float32).astype(int)[:,0])
            people_cpd.append(np.array(person['smoking_history'], dtype=np.float32).astype(int)[:,1])
        else: # never smokers
            people_ages.append(1)  # arbitrary age -- must be max_age-2 or less
            people_cpd.append(0)
        people.append(person)
    return (people,people_ages,people_cpd)

if __name__ == '__main__':
#    for year in range(1890, 2110, 10):
    year = sys.argv[2] # set year as birth_cohort passed to it in runSHG.R
    # year = 1950
    (people,people_ages,people_cpd) = parse_results(year)
    open('people_ages_mat_' + str(year) + '.txt','w').close()
    open('people_cpd_mat_' + str(year) + '.txt','w').close()
    open('people_mat_' + str(year) + '.txt','w').close()
    for item in people:
        f_handle = open('people_mat_' + str(year) + '.txt','a')
        np.savetxt(f_handle,(np.atleast_2d([item['race'], item['sex'], item['yob'], item['init_age'], item['cess_age'], item['ocd_age']])).astype(int),fmt='%i')
        f_handle.close()
    for item in people_ages:
        f_handle = open('people_ages_mat_' + str(year) + '.txt','a')
        np.savetxt(f_handle,(np.atleast_2d(item)).astype(int),fmt='%i')
        f_handle.close()
    for item in people_cpd:
        f_handle = open('people_cpd_mat_' + str(year) + '.txt','a')
        np.savetxt(f_handle,(np.atleast_2d(item)).astype(int),fmt='%i')
        f_handle.close()
