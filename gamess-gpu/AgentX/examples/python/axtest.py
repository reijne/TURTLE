#used for debug

#import methods from the python wrapper

from libpyagentx import *

#initialise

axParserStart()

#open the data documents

axDataGetUri("/home/pac83/AgentX/trunk/examples/xml/dlpoly.xml")

#open the supplementary documents

axGetUri("/home/pac83/AgentX/trunk/ontology/ontology.owl")

axGetUri("/home/pac83/AgentX/trunk/map/map.rdf")

