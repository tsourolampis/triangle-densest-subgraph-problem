//#include <bits/stdc++.h>
#include <iostream>
#include <vector>
#include <queue>          // std::priority_queue
#include <iomanip>
#include "IO.h"
#include "HASH.h"

#define MAX(a,b)) (((a)>=(b))?(a):(b))

using namespace std; 

/*
Compile: g++ TrianglePeel.cpp -o TrianglePeel -std=gnu++0x -O3
Demo: ./TrianglePeel.exe < toy.txt > toy.log
Charalampos E. Tsourakakis, babis@seas.harvard.edu 
*/

const int MAXV = 10000000;
const int MAXE = 100000000;
const int QUERYBUFFER = 6000000;

int  MAXDEG=0; 
// we shall assume that all vertices are numbered  from 1 to n 
 

int E, V;
int eList[MAXE+1][2]; 
int degrees[MAXV+1];  
int NumTriangles=0;

void SimplifyGraph() 
{
  int E1 = 1;
  int multipleEdges=0;
  HASH::init();
  for(int i = 1; i <= E; ++i) {
    //printf("Testing edge (%d,%d)\n", eList[i][0], eList[i][1]);
    if(1 <= eList[i][0] && eList[i][0] <= V &&
         1 <= eList[i][1] && eList[i][1] <= V &&
         eList[i][0] != eList[i][1] &&
         HASH::insert(eList[i][0], eList[i][1])) {
      eList[E1][0] = eList[i][0];
      eList[E1][1] = eList[i][1];
	  degrees[eList[i][0]]++;
	  degrees[eList[i][1]]++; 
      E1++;
    }
	else
		multipleEdges++;  
  }
  E = E1-1;
  
  cout<<"Number of edges in simple graph:"<<E<<endl; 
  cout<<"Number of multiple edges and self-loops that were removed:"<<multipleEdges<<endl; 

}



void GraphIn() {
  int u, v;
  IO::read(V);
  IO::read(E);
  printf("Number of vertices and edges (%d,%d)\n",V,E);
  for(int i = 1; i <= E; ++i) 
  {
    IO::read(u);
    IO::read(v);
	if(v>u)
	{
		eList[i][0] = u;
		eList[i][1] = v;
	}
    if(u>v)	
	{
		eList[i][0] = v;
		eList[i][1] = u;
	}
  }
}

void PrintEdgeList()
{
   for(int i=1; i<=E; i++) 
   {
        printf("(%d,%d)\n",eList[i][0],eList[i][1]); 
   }
}

void PrintDegrees()
{ 
   printf("***************************\n");
   printf("Vertex id \t Degree\n"); 
   for(int i = 1; i <= V; i++)
      printf("%d \t %d\n",i, degrees[i]); 
   printf("***************************\n");
}



int MaximumDegree()
{  
   for(int i = 1; i <= V; i++)
    if(MAXDEG < degrees[i] )
    	MAXDEG=degrees[i];
}


void graphinStdio() {
  scanf("%d%d", &V, &E);
  for(int i = 0; i < E; ++i) {
    scanf("%d%d", &eList[i][0], &eList[i][1]);
  }
}

void ElapsedTime(double startTime)
{ 
   printf("Elapsed time to read input: %f\n", (clock()-startTime)/CLOCKS_PER_SEC );
}

void ElapsedTime(double startTime, char* s)
{ 
   printf("Elapsed time for %s: %f\n", s, (clock()-startTime)/CLOCKS_PER_SEC );
}
///////////////adjacency list
int eStart[MAXV];
int mynext[MAXE];

const int NOEDGE = -1;

void BuildGraph() {
  for(int i = 1; i <= V; ++i) {
    eStart[i] = NOEDGE;
  }
  for(int i = 1; i <= E; ++i) {
    mynext[i] = eStart[eList[i][0]];
    eStart[eList[i][0]] = i;
  }
}

int triangles[MAXV+1]; 


void PrintTriangles()
{  
   printf("***************************\n");
   printf("Vertex id \t Triangles\n"); 
   for(int i = 1; i <= V; i++)
      printf("%d \t %d\n",i, triangles[i]); 
   printf("***************************\n");
}


void CountTrianglesNaively()
{ 
    MaximumDegree();
    int adjacent[MAXDEG];
	for(int i = 1; i <= V; ++i) 
	{
		int t = 0; 
		for( int j = eStart[i]; j!= NOEDGE; j = mynext[j] ) 
		{
			adjacent[t++] = eList[j][1];
		}
		for(int j1 = 0; j1 < t; ++j1) 
		{
			for(int j2 = 0; j2 < j1; ++j2) 
			{
				if(HASH::find(adjacent[j1], adjacent[j2]) )
				{
					triangles[i]++; 
					triangles[adjacent[j1]]++;
					triangles[adjacent[j2]]++;
                    NumTriangles += 1; 					
				}
			}
        }	
	}
	printf("Total number of triangles %d\n",NumTriangles); 
} 





int MAXTRI=-1; 

void MaximumTriangles()
{  
   for(int i = 1; i <= V; i++)
    if( MAXTRI < triangles[i] )
 	      MAXTRI=triangles[i];
}

 
vector<int> AdjList[MAXV];
 
void BuildAdjList()
{ 
  for(int i = 1; i<=E; i++)
  {
       AdjList[ eList[i][0] ].push_back( eList[i][1] ); 
	   AdjList[ eList[i][1] ].push_back( eList[i][0] );   
  }
}


int permutation[MAXV+1]; 

double  EDGEDENSITY = -1; 
int CharikarSize=-1; 
double CharikarFe = 0.0; 

 



void PrintAdjList()
{
	for(int j = 1; j<=V; j++)
	{
	    cout<< "Neighbors of "<<j<<" [";
		for (int i=0; i<AdjList[j].size();i++)
		{
		cout << AdjList[j][i]<< " ";
		}
        cout<< "]" <<endl;
		
	}
}



 


double  TRIANGLEDENSITY = -1; 
int TrianglePeelSize=-1; 
double TrianglePeelSizeFraction = -1.0;
double TrianglePeelFe = -1.0; 
double  TrianglePeelEdgeden = -1.0; 
double TrianglePeelTe = -1.0;
double TriangleOptVals[MAXV+1]; 

int OPTIND; 


priority_queue< pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>> > q;

void PQPeelTriangles()
{

    TRIANGLEDENSITY= 3.0*NumTriangles/V; 
	TrianglePeelSize=V;
	TrianglePeelSizeFraction=1.0;
	TrianglePeelFe = 2*E/(V*(V-1));
	TrianglePeelTe = 6*NumTriangles/(V*(V-1)*(V-2));
	TrianglePeelEdgeden = 2*E/V;
	
	double numedges=(double)E;
	double numoftriangles=NumTriangles;
	double numvertices=(double)V;	
	
	
	
	for (int i = 1; i <= V; ++i) 
        q.push(make_pair(triangles[i], i));

	int c = q.top().second;
	int counter = 0; 
	while (!q.empty()) 
	{
		int c = q.top().second;
        q.pop();
         if (triangles[c] < 0 )// == -1) 
		{ 
            continue;
        }
		else
		{ 
			TriangleOptVals[counter] = 3.0*numoftriangles/numvertices ;
			permutation[counter++]=c; 
		}
		
        numedges -= degrees[c]; // number of edges goes down by degree of c
        --numvertices; // one vertex less now 
		numoftriangles -= triangles[c];
        triangles[c] = -1; //no need to look at it again 
		if (numvertices > 0) 
		{
			if( TRIANGLEDENSITY <= 3.0*numoftriangles/numvertices  ) 
			{
				TrianglePeelEdgeden = 2.0*numedges/numvertices;
				TRIANGLEDENSITY = 3.0*numoftriangles/numvertices ;
				TrianglePeelSize = numvertices;
				TrianglePeelSizeFraction = numvertices/V;
				TrianglePeelFe = 2*numedges/(numvertices*(numvertices-1));
				TrianglePeelTe = 6*numoftriangles/(numvertices*(numvertices-1)*(numvertices-2));
				OPTIND = counter; 
				//subgraph.clear();
			}
        } 
		for (int i=0; i < AdjList[c].size();i++)
		{
			int w = AdjList[c][i];
			degrees[w]--;
			int Tw = triangles[w];
			if( Tw >= 0 && degrees[c] >=2 ) 
			{
			for(int j = 0; j < AdjList[w].size(); j++)
			{	
  				    int candidate = AdjList[w][j];					
					if( triangles[candidate] > 0 && candidate !=c)
					{
					    
						if( degrees[c] <= degrees[candidate] )
						{
							for (int rc=0; rc < AdjList[c].size(); rc++)
							{
							  if( AdjList[c][rc] == candidate)
									triangles[w]--; 		  								 						
							}
						}
						else
						{
							for (int rc=0; rc < AdjList[candidate].size(); rc++)
								if( AdjList[candidate][rc] == c)									
											triangles[w]--; 		
						}
						
					}				
					
				}
				
				q.push(make_pair(triangles[w], w));				
			}	
			 
		}
		if( counter == V )
		    break;
	}

}

void PrintPermutationStatistics()
{  
   int breakline = 0; 
   cout<<"Peeling permutation order"<<endl; 
   cout<<"[";
   for(int i=0; i < V; i++) 
   {
		cout<<" "<< permutation[i];
		if( ++breakline > 25 ) // entries per line 
		{
		        cout<<endl;
				breakline = 0;
		}
   }
   cout<<"]"<<endl;
  
   
   breakline = 0;
   cout<<"Optimal triangle density values"<<endl;
   cout<<"[";
   for(int i=0; i < V; i++) 
   {
		cout<<" "<< TriangleOptVals[i];
		if( ++breakline > 25 ) // entries per line 
		{
		        cout<<endl;
				breakline = 0;
		}
   }
   cout<<"]"<<endl;
   
}



void printOptSubgraph()
{ 
   cout<<"[";
   for(int i=OPTIND; i < V; i++)  
   {
		cout<<" "<< permutation[i];
   }
   cout<<"]"<<endl;

}





int main(int argc, char **argv)
{
  double startTime = clock();
  GraphIn();
  ElapsedTime(startTime);
  //PrintEdgeList();
  SimplifyGraph(); 
  //PrintDegrees();
  //MaximumDegree();
  BuildGraph();
  BuildAdjList();
  //PrintAdjList();
  
  /* Count first */ 
  cout<<"**************** Triangle  Counting *****************"<<endl;
  startTime = clock();
  CountTrianglesNaively();
  ElapsedTime(startTime,"Triangle Counting"); 
  
  cout<<"******************************************************"<<endl;
  startTime = clock();
  PQPeelTriangles();
  ElapsedTime(startTime,"Triangle Peeling");
  cout<<"************** Triangle  PEELING ***************"<<endl;
  cout<<"Triangle Peeling's Results"<<endl; 
  cout<<"Optimal avg. triangle density:"  <<TRIANGLEDENSITY<<endl;
  cout<<"Size:"  <<setprecision(3) << fixed <<TrianglePeelSize<<endl;
  cout<<"Size/|V|:"  <<setprecision(3) << fixed <<TrianglePeelSizeFraction<<endl;
  cout<<"fe:"  <<setprecision(3) << fixed << TrianglePeelFe<<endl;
  cout<<"*******************************************"<<endl;
  bool subgraphflag = false; //false 
  if( subgraphflag ) 
  {
      cout<<"********** Optimal Subgraph *********************"<<endl;
      printOptSubgraph();
	  cout<<"*************************************************"<<endl;
  }

  bool statsflag = false; //false 
  if( statsflag ) 
  {
     
	  cout<<"********** Permutation Statistics *********************"<<endl;
	  PrintPermutationStatistics();
	  cout<<"*************************************************"<<endl;
  }
  
  

   
  /* this works too but is slower, kept it for sake of completeness
  startTime = clock();
  PeelTriangles();
  cout<<"************** Triangle  PEELING ***************"<<endl;
  ElapsedTime(startTime,"Triangle Peeling"); 
  cout<<"Triangle Peeling's Results"<<endl; 
  cout<<"Avg. Triangle density:"  <<TRIANGLEDENSITY<<endl;
  cout<<"Size:"  <<TrianglePeelSize<<endl;
  cout<<"fe:"  <<TrianglePeelFe<<endl;
  cout<<"ft:"  <<TrianglePeelTe<<endl;
  cout<<"*******************************************"<<endl;
  ElapsedTime(startTime,"Triangle peeling"); 
  */
 
  return 0;
}
 
