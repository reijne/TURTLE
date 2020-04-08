
/*==========================================================================
    FILE          : Exp_eval.CPP
  
    AUTHOR        : Ashish & Parth
    
    CONTECT       : parth2651@yahoo.co.in

	 CREATED       : 15 March 2005
    
	 LAST MODIFIED : 15 March 2005

    This Program Parses an algebraic expression.Converts it to Postfix
	 and then carries out evaluation of the postfix expression.

==========================================================================*/

/* Modified by P Couch 04/05/05 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<ctype.h>
#include<string.h>

#ifdef HAVE_CONFIG_H
#include<config.h>
#endif

static char* eq=NULL;
static char stack[20],temp[32],a;
static int topv=0,i=0,j=0,k,l,m,n=0,top=0,p=0,braces=0;
static float stackv[20];
static char post[50];

static struct
  {
    char name[32];
    float value;
    int isvar;
  }symb_table[20];

void push(char t)
  {
    stack[top]=t;
	 top++;
  }

char pop()
  {
    top--;
    return stack[top];
  }

int emptystack()
  {
    if(top==0)
      return 1;
	 else
      return 0;
  }

int precedence(char temp)
  {
    if(temp=='('||temp==')')
      return 3;
    if(temp=='^')
      return 2;
    if(temp=='*'||temp=='/'||temp=='%')
      return 1;
    else
		return 0;
  }

int isoperator(char t)
  {
    if(t=='+'||t=='-'||t=='*'||t=='/'||t=='%'||t=='^')
      return 1;
    else
      return 0;
  }

int isoperand(char t)
  {
	 if(!isoperator(t)&&t!=')'&&t!='('&&t!='\0')
      return 1;
    else
      return 0;
  }

void addtemp(char t)
  {
    if(m==0)
      j=0;
    temp[j]=t;
      j=j+1;
    temp[j]='\0';
  }

int isreal(char temp[])
  {
    int l,j=0,count=0;
    l=strlen(temp);
    for(j=0;j<l;j++)
      {
	if(temp[j]=='.')
	  count++;
	if(!((temp[j]>=48&&temp[j]<=57)||temp[j]=='.'))
	  return 0;
      }
  if(count>1)
	 {
		printf("ERROR : there can not be two '.' in one no.");
		return(0);
	 }
  return 1;
  }

int validoperand(char temp[])
  {
    int a=0;
    char t[32];
	 k=0;
	 l=strlen(temp);
	 while(temp[k]==' ')
		k++;
	 while(temp[l-1]==' ')
		l--;
	 for(a=k;a<l;a++)
		t[a-k]=temp[a];
	 t[a-k]='\0';
	 strcpy(temp,t);
	 l=strlen(temp);
	 if(!isreal(temp)&&isdigit(temp[0]))
		{
	printf("ERROR : first ch of operand must be a character ");
	return(0);
		}
	 for(k=0;k<l;k++)
		{
	if(!isalnum(temp[k])&&temp[k]!='.'&&temp[k]!='_')
	  {
		 return 0;
	  }
		}
	 return 1;
  }

void addsymboltable(char temp[])
  {
	 strcpy(symb_table[n].name,temp);
	 if(isreal(temp))
		{
	int j=0,l,k=0;
	float t=0;
	l=strlen(temp);
	for(j=0;j<l;j++)
	  {
		 if(temp[j]=='.')
			{
		k=1;
		j++;
			}
		 if(k==0)
			t=t*10+(temp[j]-48);
		 else
			{
		float t1=pow(10,k*(0-1));
		t=t+(temp[j]-48)*t1;
		k++;
			}
	  }
	symb_table[n].value=t;
	symb_table[n].isvar=0;
		}
	 else
		{
	symb_table[n].value=0;
	symb_table[n].isvar=1;
		}
	 n++;
  }

int sementicrules(int p)
  {
	 switch (p)
		{
		case 1:    /*for operator*/
	  if(isoperator(eq[i+1])||eq[i+1]==')'||eq[i-1]=='(')
		 return 0;
		 break;
		case 2:     /*for ')'*/
		  if((eq[i+1]=='(')||((isoperand(eq[i+1]))&&(eq[i+1]!=' '))||(isoperator(eq[i-1]))||(eq[i-1]=='('))
		 return 0;
		 break;
		case 3:    /*for operand*/
	  if(eq[i+1]=='(')
		 return 0;
		 break;
		}
	 return 1;
  }

int parse_infix(char eq[])
  {
	 i=0;
	 while(eq[i]==' ')
		i++;
	 if(isoperator(eq[i])||eq[i]==')')
		{
	printf("\nERROR : first ch can not be operator or ')' ");
	return(0);
		}
	 while(eq[l-1]==' ')
		l--;
	 eq[l]='\0';
	 if(isoperator(eq[l-1]))
		{
	printf("ERROr : last ch can not be an operator");
	return(0);
		}
	 while(eq[i]!='\0')
		{
	a=eq[i];
	if(isoperator(a))
	  {
		 if(!sementicrules(1))
			{
		printf("\nERROR : invalid use of operator  ");
		return(0);
			}
		 while(eq[i+1]==' ')
			i++;
		 if(isoperator(eq[i+1])||eq[i+1]==')')
			{
		printf("\nERROR : invalid use of operator  ");
		return(0);
			}
	  }
	else if(a=='(')
	  {
		 braces=braces+1;
		 while(eq[i+1]==' ')
			i++;
	  }
	else if(a==')')
	  {
		 if(braces>0)
			{
		braces=braces-1;
		if(!sementicrules(2))
		  {
			 printf("\nERROR : invalid use of ')' ");
			 return(0);
		  }
			}
		 else
			{
		printf("\nERROR : no of opening & closing braces should be equal ");
		return(0);
			}
		 while(eq[i+1]==' ')
			i++;
	  }
	else /*if(isoperand(a))*/
	  {
		 m=0;
		 while(isoperand(a))
			{
		if(!sementicrules(3))
		  {
			 printf("\nERROR : invalid use of operand ");
			 return(0);
		  }
		addtemp(a);
		m=1;
		i++;
		a=eq[i];
			}
		 if(validoperand(temp))
			addsymboltable(temp);
		 else
			{
		printf("\nERROR : invalid operand ");
		return(0);
			}
		 if(m)
			i--;
	  }
	i++;
		}
	 if(braces!=0)
		{
	printf("\nERROR : no. of opening & closing braces should be equal");
	return(0);
		}
	 return(0);
  }

void addpost(char p1)
  {
	 if(p1!=' ')
		{
	post[p]=p1;
	p=p+1;
	post[p]='\0';
		}
  }

void in_post(char eq[])
  {
    char t1;
	 i=0;
	 a=eq[i];
	 push('(');
	 eq[l++]=')';
	 eq[l]='\0';
	 while(a!='\0')
		{
	if(a=='(')
	  push(a);
	if(isoperand(a))
	  addpost(a);
	if(isoperator(a))
	  {
		 while(!emptystack()&&(precedence(stack[top-1])>=precedence(a)))
			{
		if(stack[top-1]=='(')
		  break;
		t1=pop();
		if(t1!='(')
		  addpost(t1);
			}
		 push(a);
	  }
	if(a==')')
	  {
		 char t1=pop();
		 while(t1!='(')
			{
		addpost(t1);
		t1=pop();
			}
		 push(t1);
	  }
	i++;
	a=eq[i];
		}
	 while(top>0)
		{
	char t1=pop();
	if(t1!='(')
	addpost(t1);
		}
  }

void pushv(float tt)
  {
	 stackv[topv]=tt;
	 topv++;
  }

float popv()
  {
	 topv--;
	 return stackv[topv];
  }

float operation(float v1,float v2,char t)
  {
	 if(t=='*')
		return (v1*v2);
	 if(t=='/')
		return (v1/v2);
	 if(t=='+')
		return (v1+v2);
	 if(t=='-')
		return (v1-v2);
	 /*	 if(t=='%')
		 return (int(v1)%(int(v2)));*/
	 if(t=='^')
		return pow(v1,v2);
	 return 0;
  }

float evaluate(char post[])
  {
	 int i=0,j,k,l1,p=0;
	 float v1,v2;
	 for(k=0;k<l;k++)
		{
	if(isoperator(post[k]))
	  {
		 v2=popv();
		 v1=popv();
		 v1=operation(v1,v2,post[k]);
		 pushv(v1);
	  }
	else
	  {
		 l1=strlen(symb_table[i].name);
		 if(symb_table[i].isvar==1)
			{
		p=0;
		for(j=0;j<i;j++)
		  {
			 if(strcmp(symb_table[i].name,symb_table[j].name)==0&&symb_table[j].value>0)
				{
			symb_table[i].value=symb_table[j].value;
			p=1;
			break;
				}
		  }
		if(p==0)
		  {
			 printf("Enter value of %s : ",symb_table[i].name);
			 scanf("%f",&symb_table[i].value);
		  }
			}
		 pushv(symb_table[i].value);
		 k=k+l1-1;
		 i++;
	  }
		}
	 return popv();
  }

float evaluateExpression(char* expression){

  eq=expression;
  
  l=strlen(eq);
  parse_infix(eq);
  l=strlen(eq);
  in_post(eq);
  l=strlen(post);
  
  /* reset statics */

  eq=NULL;
  topv=0;
  i=0;
  j=0;
  n=0;
  top=0;
  p=0;
  braces=0;

  return evaluate(post);
  
}
