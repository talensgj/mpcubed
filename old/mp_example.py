#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import multiprocessing
import os
import random
import time

class simulation(multiprocessing.Process):
    def __init__(self, name):
        # must call this before anything else
        multiprocessing.Process.__init__(self)

        # then any other initialization
        self.name = name
        self.number = 0.0
        sys.stdout.write('[%s] created: %f\n' % (self.name, self.number))

    def run(self):
        sys.stdout.write('[%s] running ...  process id: %s\n' 
                         % (self.name, os.getpid()))
        time.sleep(10)
        self.number = random.uniform(0.0, 10.0)
        sys.stdout.write('[%s] completed: %f\n' % (self.name, self.number))

sim_list = []
sim_list.append(simulation('foo'))
sim_list.append(simulation('bar'))

for sim in sim_list:
    sim.start()
