
/*simple c example*/

#include<stdio.h>
#include<agentx.h>

int main(){
  
  int i=0,j=0,k=0,noBasisSetAssignment=0,noMolecule=0,noAtom=0,noAtomicBasisSet=0,noBasisGroup=0;
  int noMolecularOrbital=0;
  int noprop=0;
 
  char **results = NULL;

  /* initialise parser */

  if ( axParserStart() )
    {
      fprintf(stderr, "Error initialising AgentX\n");
    }

  /* parse AgentX control file (AgentX.ini) */

  if ( axControl("agentx.ini",(NULL)) )
    {
      fprintf(stderr, "Error parsing agentx.ini\n");
    }

noprop = axPath( "Molecule.Atom[*].xCoordinate,yCoordinate,zCoordinate", &results );  

  noMolecule=axSelect("Molecule");
  
  for(j=0;j<noMolecule;j++){
    
    noprop=axSelect("identifier");
    if(noprop>0){
      printf("\nserialisation of Molecule :%s",axValue());
      axDeselect();      
    }
    
    noAtom=axSelect("Atom");
    
    for(k=0;k<noAtom;k++){
      
      noprop=axSelect("identifier");
      if(noprop>0){
	printf("\nserialisation of Atom :%s\n",axValue());
	axDeselect();
      }
      
      noprop=axSelect("elementType");
      if(noprop>0){
	printf("element = %s\n",axValue());
	axDeselect();
      }
      
      noprop=axSelect("xCoordinate");
      if(noprop>0){
	printf("x coordinate = %s\n",axValue());
	axDeselect();
      }
      
      noprop=axSelect("yCoordinate");
      if(noprop>0){
	printf("y coordinate = %s\n",axValue());
	axDeselect();
      }
      
      noprop=axSelect("zCoordinate");
      if(noprop>0){
	printf("z coordinate = %s\n",axValue());
	axDeselect();
      }
      
      axSelectNext();
      
    }
    axDeselect();
    axSelectNext();
  }
  
  axDeselect();
  
  axParserFinish();
  
  return (0);
  
}    
