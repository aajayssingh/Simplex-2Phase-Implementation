/*
two phase simplex
Assuming simplex tableau wlog of form
	x1	x2	s1	s2	a1 a2rhs                    //note order of vars, slacks and auxiliary vars
z                                               // note all z have all vars shifted left
s1
s2

slacks are x3 and x4

> since input is in initial tableu from standard form all rhs is supposed to be non negative.
> since input is initial tableau from standard form so we are supposed to solve for maximization problem.
> constrainsts are assumed to be  <= form.

phase 1    //http://nsdl.niscair.res.in/jspui/bitstream/123456789/826/1/CHAPTER-2%20LP.pdf

>Max prob
    with most negative as entering column.
> find identity matrix for initial basis vars.
> and modify objective is in input
>another assumption is that # unknowns = # constraints = m-1

It should be noted that although Case z* = 0 implies that all the artificial variables are equal to zero, this does not mean that they are all out of the basis.

So it is necessary to consider it in more detail, namey:

Case 2.1: All the artificial variables are non-basic.

Case 2.2: some of the artificial variables are in the basis




*/


#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <iomanip>

using namespace std;
double **coeffMatrix = 0;
double **coeffMatrixPhase1 = 0;
int m = 0,n = 0, a = 0, ukn = 0;
int searchMinz();

//calculates leaving basis vector index on the basis of min-ratio
//sanity checks to ensure min ratio is positive only.
int calMinRatio_LeavingBasisVector(int enteringVectorIndex);
void displayCoeffiecientMatrix();

/*
desp :  prints the current basis variable values
author : Ajay Singh
*/
void printSolutionVector();

/*
desp :  update the basis
author : Ajay Singh
*/
void updateSolutionvector(int leavingIndex, int enteringIndex);

void displayPhase1Basis();

/*
desp :  display tableau
author : Ajay Singh
*/
void displayCoeffiecientMatrixPhase1();
int findIdentityMatrixPhase1();

/*
desp :  execute the phase 1 of simplex algorithm
author : Ajay Singh
*/
void simplexPhase1();

/*
desp :  after entering and leaving of variables operation update the basis vectorfor phase1
author : Ajay Singh
*/
void updateSolutionvectorPhase1(int leavingIndex, int enteringIndex);
int calMinRatio_LeavingBasisVectorPhase1(int enteringVectorIndex);
int searchMinzPhase1();

/*
desp :  execute the phase 2 of simplex algorithm
author : Ajay Singh
*/
void simplexPhase2();

/*
desp :  after entering and leaving of variables operation update the basis vectorfor phase2
author : Ajay Singh
*/
void updateSolutionvectorPhase2(int leavingIndex, int enteringIndex);

/*
desp :  finds if the colIndex is on amongst the artificial variables
author : Ajay Singh
*/
bool artColinBasis(int colIndex);//col index starts from 1

bool phase1flag1 = false;
bool basisContainsArtificial = false;

vector<pair<int, double>> basisVector ;
vector<pair<int, double>> basisVectorPhase1 ;//int = x_i th var from 1 but matrix has zero so take care

int main()
{
    int i = 0,j = 0;

    ifstream input_file;
    string delimiter1 = " ", str_item, inputLine;
    size_t pos = 0;
    string::size_type sz;

    input_file.open("inp-params.txt", ifstream::in);

    cout << "maximization standard problem " <<endl;
    cout << "enter m(#rows in std input) : "; cin>>m;
    cout << "enter n(#cols in std input) : "; cin>>n;
    cout << "#artificial vars(placed at last cols) :"; cin>>a;
    cout << "#unknowns ie x vars : "; cin>>ukn;

    // dynamically allocated coefficient matrix for simplex tableau
    //solutionVector = new double* [n];
    coeffMatrix = new double* [m];
    coeffMatrixPhase1 = new double* [m];

    for(i = 0; i < m; i++)
    {
        coeffMatrix[i] = new double [n];
        coeffMatrixPhase1[i] = new double [n];
    }


    for(int i = 0 ; i < m ; i++ )
    {
        cout<<endl;
        for( int j = 0 ; j < n ; j++ )
        {
            coeffMatrix[i][j] = 0;
            coeffMatrixPhase1[i][j] = 0;
        }
    }


    if (input_file.is_open())
    {
        for( int i = 0 ; i < m ; i++ )
        {
            cout << endl;
            if(getline(input_file, inputLine))
            {
                for( int j = 0 ; j < n ; j++ )
                {
                    pos = inputLine.find(delimiter1);
                    str_item = inputLine.substr(0, pos);
                    //cout<<str_item<< endl;
                    coeffMatrix[i][j] = stod(str_item, &sz);
                    inputLine.erase(0, pos + delimiter1.length());
                }
            }

        }

    }
    //copy  coeff matrix to phase1 matrix
    for( int i = 0 ; i < m ; i++ )
    {
        for( int j = 0 ; j < n ; j++ )
        {
            if((0 == i)&&(j < (n-a-1)))
            {
                //all vars except artificial are 0
            }
            else{
            coeffMatrixPhase1[i][j] = coeffMatrix[i][j];
            }
        }
    }

    displayCoeffiecientMatrix();
    findIdentityMatrixPhase1();

    cout <<"modified tableu :" <<endl ;
    displayCoeffiecientMatrixPhase1();
    displayPhase1Basis();

    simplexPhase1();

    //postphase1 procesing and case handling :
    if(phase1flag1)
    {
        if(coeffMatrixPhase1[0][n-1] == 0)//Z* = 0
        {
            //check if basis has artificial var or not
            for(int i = 0; i < basisVectorPhase1.size(); i++)
            {
                if((basisVectorPhase1[i].first - 1 >= (n-a-1) ))//xi starts from 1
                {
                    cout <<"basis contains artificial x"<<basisVectorPhase1[i].first<<" = "<<basisVectorPhase1[i].second<<endl;
                    basisContainsArtificial = true;
                }
                //
            }
            if(basisContainsArtificial)
            {
                cout<<"phase 1 feasible but contains artificial vars"<<endl;
            }
            else{
            cout<<"phase 1 feasible No artificial vars"<<endl;
            }
            cout<<"move to phase 2";
            simplexPhase2();
        }
        else{
        //If all zj – cj ≥ 0 and at least one artificial variable occurs in the optimum basis and hence max z* < 0,
        //then LPP has infeasible solution.
        cout<<"phase 1 Not feasible!!!!"<<endl;
        }
    }

     //free dynamically allocated coefficient matrix for tableau
    for( int i = 0 ; i < m ; i++ )
        delete[] coeffMatrix[i] ;

    delete [] coeffMatrix;
    return 0;
}

//calculates min reduced cost, to decide the index of entering column to basis, returns index of min reduced cost column
//user shall check if its positive or negative to decide optimal solution is reached or not.
int searchMinz()
{
    int min_index = 0;
    int min_val = coeffMatrix[0][0];
    int j;
    for(j = 0; j < n -1; j++)
    {
        if(coeffMatrix[0][j] < 0)
        {
            break;//at least one negative z element
        }
    }

    if(j == n - 1)
    {
         cout <<"optimal solution reached !"<<endl;
         return -1;
    }

    //else there is atleast one negative ele

    for(int i = 1; i < n -1; i++)
    {
        if(min_val > coeffMatrix[0][i])
        {
            min_val = coeffMatrix[0][i];
            min_index = i;
        }
    }

    if(min_val < 0)
    {
        int i = 0;
        //check if atleast one column value is positive , if not solution is unbounded
        for(i = 1; i < m; i++)
        {
            if(coeffMatrix[i][min_index] > 0)
                break;
        }

        if(i == m)
        {
        cout <<"solution is unbounded "<<endl;//all elements in column are negative
        return -1;
        }
    }
    return min_index;

}

//calculates leaving basis vector index on the basis of min-ratio
//sanity checks to ensure min ratio is positive only.
int calMinRatio_LeavingBasisVector(int enteringVectorIndex)
{
    double min_ratio_val = 999999, new_ratio;//default large min ratio
    int min_ratio_index = -1;//default row index;

    for(int i = 1; i < m; i++)//iterate all constraint rows
    {
        if(coeffMatrix[i][enteringVectorIndex] != 0)
        {
            new_ratio = coeffMatrix[i][n-1]/coeffMatrix[i][enteringVectorIndex];
            if(new_ratio < 0)
            {
                cout<<"min ratio negative continue to cal next min ratio"<<endl;
                continue;
            }
            else if(((min_ratio_val) > new_ratio ))
            {
                min_ratio_val = new_ratio;
                min_ratio_index = i;
            }
            else{
            cout<<"current rows min ratio not negative nor is minimum "<<new_ratio<<endl;
            }

        }
    }

    if(min_ratio_val >= 0)
    {
        return min_ratio_index;//return leaving  basis row;
    }
    else
    {
        cout<<"all min ratios are negative";
        return -1;// all min ratios are negative
    }
}

//update the solution vector at every iteration
void updateSolutionvector(int leavingIndex, int enteringIndex)
{
    basisVector[leavingIndex - 1].first = enteringIndex+1;//stores x1, x2 from 0th index of vector
    basisVector[leavingIndex - 1].second = coeffMatrix[leavingIndex][n-1];
    for(int i = 0; i < basisVector.size(); i++)
    {
        if(basisVector[i].first != enteringIndex+1)
        {
            basisVector[i].second = coeffMatrix[i+1][n-1];
        }
    }

}

/*
desp :  display the current basis values
author : Ajay Singh
*/
void printSolutionVector()//non phase1
{
    for(int i = 0; i < basisVector.size(); i++)
    {
        cout <<"x"<<basisVector[i].first<<" = "<<basisVector[i].second<<endl;

    }
}

/*
desp :  display tableau
author : Ajay Singh
*/
void displayCoeffiecientMatrix()
{
    for(int i = 0 ; i < m ; i++ )
    {
        cout<<endl;
        for( int j = 0 ; j < n ; j++ )
        {
           cout<< coeffMatrix[i][j] <<"   ";
        }
    }
}

/*
desp :  find current basis vars which are in Identity matrix form
author : Ajay Singh
*/
int findIdentityMatrixPhase1()
{
//needs #unknows = m-1
int countNonzero = 0;
int Numslacks = n- a -ukn -1;
cout<<endl<< "#slacks = " <<Numslacks<<endl;
    for(int i = ukn; i < n-1 ; i++)//each col except last as it has nothing for identity chkin
    {
        for(int j = 1; j < m; j++) //evry row inc z
        {
            if(i < (ukn + Numslacks))//for slack/surplus cols
            {
                if(coeffMatrixPhase1[j][i] == 1)
                {
                    basisVectorPhase1.push_back(pair<int, double> (i + 1, coeffMatrixPhase1[j][n-1]));
                    break;
                }
            }
            else{//for artificial vars
                    if(coeffMatrixPhase1[j][i] == 1)
                    {
                        basisVectorPhase1.push_back(pair<int, double> (i + 1, coeffMatrixPhase1[j][n-1]));

                        //make tableau a proper tableau
                        for(int k = 0; k < n; k++)//each column k of rowj effect row 0
                        {
                            coeffMatrixPhase1[0][k] = coeffMatrixPhase1[0][k] - (coeffMatrixPhase1[j][k]);
                        }


                    }

            }
        }

    }

    for(int i = 0; i < basisVectorPhase1.size(); i++)
    {
        //cout <<"x"<<basisVectorPhase1[i].first<<" = "<<basisVectorPhase1[i].second<<endl;

        basisVector.push_back(pair<int, double> (basisVectorPhase1[i].first, basisVectorPhase1[i].second));
    }
    cout<<"before basis vectorPh1 :"<<endl;
    displayPhase1Basis();

    for(int i = 0; i < basisVectorPhase1.size(); i++)
    {
        //cout <<"x"<<basisVectorPhase1[i].first<<" = "<<basisVectorPhase1[i].second<<endl;
        int colIndex = basisVector[i].first - 1;
        for(int j = 1; j <m; j++)
        {
            if(coeffMatrixPhase1[j][colIndex] == 1)
            {
                basisVectorPhase1[j-1].first = basisVector[i].first;
                basisVectorPhase1[j-1].second = basisVector[i].second;
                break;
            }
        }
    }

    cout<<"after ordering basis vectorPh1 :"<<endl;
    displayPhase1Basis();


}

void displayPhase1Basis()
{
    for(int i = 0; i < basisVectorPhase1.size(); i++)
    {
        cout <<"x"<<basisVectorPhase1[i].first<<" = "<<basisVectorPhase1[i].second<<endl;
    }
}


/*
desp :  display tableau
author : Ajay Singh
*/
void displayCoeffiecientMatrixPhase1()
{
    for(int i = 0 ; i < m ; i++ )
    {
        cout<<endl;
        for( int j = 0 ; j < n ; j++ )
        {
           cout<< coeffMatrixPhase1[i][j] <<"   ";
        }
    }
    cout<<endl;
}

/*
desp :  execute the phase 1 of simplex algorithm
author : Ajay Singh
*/
void simplexPhase1()
{
    while(1)
    {
        int entering_index =  searchMinzPhase1();
        if((entering_index == -1)||(coeffMatrixPhase1[0][entering_index] > 0))
        {
            cout<< "stopping simplex iterations" << endl;
            break;//optimal or unbounded solution
        }
        int leaving_index = calMinRatio_LeavingBasisVectorPhase1(entering_index);

        if(leaving_index == -1)
        {
            cout<< "stopping simplex iterations" << endl;
            break;
        }

        double pivot_element = coeffMatrixPhase1[leaving_index][entering_index];

        //make pivot element unity along with row operation by dividing whole row by pivot element
        if(pivot_element != 0){
            for(int j = 0; j < n; j++)
            {
                coeffMatrixPhase1[leaving_index][j] /= pivot_element;
            }
        }
        else{
            cout<< "pivot element is 0"<< endl;
        }

        pivot_element = coeffMatrixPhase1[leaving_index][entering_index];//resetting the pivot element to the value after the roe operation ie 1
        cout<<"pivot = "<<pivot_element<<endl;
        cout<<"leaving row = "<<leaving_index<<" "<<"entering column = "<<entering_index<<endl;
        //make all elements of pivot column 0 expect pivot ele, Ri = Ri - Rp*Ri
        for( int i = 0 ; i < m ; i++ )
        {
            double multiply_factor = coeffMatrixPhase1[i][entering_index];
            if(leaving_index != i)
            {
                for(int j = 0; j < n; j++)
                {
                    coeffMatrixPhase1[i][j] = coeffMatrixPhase1[i][j] - ((multiply_factor) * (coeffMatrixPhase1[leaving_index][j]));
                    if(fabs(coeffMatrixPhase1[i][j]) < 0.000001)
                    {
                        cout<< "setting 0  "<<coeffMatrixPhase1[i][j]<<endl;
                        coeffMatrixPhase1[i][j] = 0;

                    }

                }
            }
        }

        updateSolutionvectorPhase1(leaving_index, entering_index);
        cout <<"simplex tableau at current iteration : "<<endl;

        displayCoeffiecientMatrixPhase1();
        cout<<endl;
        displayPhase1Basis();
        cout<<endl<<"********************************************************************************"<<endl;
    }


}

//update the solution vector at every iteration for phase 1
/*
desp :  after entering and leaving of variables operation update the basis vectorfor phase1
author : Ajay Singh
*/
void updateSolutionvectorPhase1(int leavingIndex, int enteringIndex)
{


    basisVectorPhase1[leavingIndex - 1].first = enteringIndex+1;//stores x1, x2 from 0th index of vector
    basisVectorPhase1[leavingIndex - 1].second = coeffMatrixPhase1[leavingIndex][n-1];
    for(int i = 0; i < basisVectorPhase1.size(); i++)
    {
        if(basisVectorPhase1[i].first != enteringIndex+1)
        {
            basisVectorPhase1[i].second = coeffMatrixPhase1[i+1][n-1];
        }
    }
}

/*
desp :  after entering and leaving of variables operation update the basis vectorfor phase2
author : Ajay Singh
*/
void updateSolutionvectorPhase2(int leavingIndex, int enteringIndex)
{
///basis vector of phase 1 is used but same coeffmatrix is used in phase2

    basisVectorPhase1[leavingIndex - 1].first = enteringIndex+1;//stores x1, x2 from 0th index of vector
    basisVectorPhase1[leavingIndex - 1].second = coeffMatrix[leavingIndex][n-1];
    for(int i = 0; i < basisVectorPhase1.size(); i++)
    {
        if(basisVectorPhase1[i].first != enteringIndex+1)
        {
            basisVectorPhase1[i].second = coeffMatrix[i+1][n-1];
        }
    }
}

int searchMinzPhase1()
{
    int min_index = 0;
    int min_val = coeffMatrixPhase1[0][0];
    int j;
    for(j = 0; j < n-1; j++)
    {
        if(coeffMatrixPhase1[0][j] < 0)
        {
            break;//at least one negative z element
        }
    }

    if(j == n-1)
    {
         cout <<"Phase 1 optimal solution reached !"<<endl;
         phase1flag1 = true;
         return -1;
    }

    //else there is atleast one negative ele


    for(int i = 1; i < n-1; i++)//n-1 so as to not to compare with value of z as last col is z
    {
        if(min_val > coeffMatrixPhase1[0][i])
        {
            min_val = coeffMatrixPhase1[0][i];
            min_index = i;
        }
    }

    if(min_val < 0)
    {
        int i = 0;
        //check if atleast one column value is positive , if not solution is unbounded
        for(i = 1; i < m; i++)
        {
            if(coeffMatrixPhase1[i][min_index] > 0)
                break;
        }

        if(i == m)
        {
        cout <<"Phase 1 solution is unbounded "<<endl;//all elements in column are negative
        return -1;
        }
    }
    return min_index;

}

//calculates leaving basis vector index on the basis of min-ratio
//sanity checks to ensure min ratio is positive only.
int calMinRatio_LeavingBasisVectorPhase1(int enteringVectorIndex)
{
    double min_ratio_val = 999999, new_ratio;//default large min ratio
    int min_ratio_index = -1;//default row index;

    for(int i = 1; i < m; i++)//iterate all constraint rows
    {
        if(coeffMatrixPhase1[i][enteringVectorIndex] != 0)
        {
            new_ratio = coeffMatrixPhase1[i][n-1]/coeffMatrixPhase1[i][enteringVectorIndex];
            if(new_ratio < 0)
            {
                continue;
            }
            else if(((min_ratio_val) > new_ratio ))
            {
                min_ratio_val = new_ratio;
                min_ratio_index = i;
            }
            else{
            cout<<"current rows min ratio not negative nor is minimum "<<new_ratio<<endl;
            }

        }
    }

    if(min_ratio_val >= 0)
    {
        return min_ratio_index;//return leaving  basis row;
    }
    else
    {
        cout<<"all min ratios are negative";
        return -1;// all min ratios are negative
    }

}

/*
desp :  execute the phase 2 of simplex algorithm
author : Ajay Singh
*/
void simplexPhase2()
{
if(!basisContainsArtificial)
{
    //copy orignal z to coeffmatrix
    for(int j = 0 ;j < n; j++)
    {
        coeffMatrixPhase1[0][j] = coeffMatrix[0][j];
    }
    //remove artificial vars columns
    //copy n-1th col to first artificial var col and reduce n
    //first artificial var col n-a-1
    for(int i = 0; i < m; i++)
    {
        coeffMatrixPhase1[i][n-a-1] = coeffMatrixPhase1[i][n-1];
    }
     //now set new n
     n = n-a;

         //to work with coeffmatrix
        for( int i = 0 ; i < m ; i++ )
        {
            for( int j = 0 ; j < n ; j++ )
            {

                coeffMatrix[i][j] = coeffMatrixPhase1[i][j] ;

            }
        }

    displayCoeffiecientMatrix();
    cout<<endl;
    displayPhase1Basis();//using basis vector of phase 1

}
else//z* = 0 but contains a artificial var with val 0
{
    cout<< "phase 2 with artificial vars"<<endl;
    //copy orignal z to coeffmatrix
    for(int j = 0 ;j < n; j++)
    {
        coeffMatrixPhase1[0][j] = coeffMatrix[0][j];
    }

    //remove artificial vars columns which is not in basis
    //set whole column to 0
    for(int j = (n-a-1) ;j < n-1; j++)
    {
        if(!artColinBasis(j+1))//var indexed from 1 is sent
        {
        //set all rows of this col = 0
            for(int i = 0 ;i < m; i++)
            {
                coeffMatrixPhase1[i][j] = 0;
            }
        }

    }

    //to work with coeffmatrix
    for( int i = 0 ; i < m ; i++ )
    {
        for( int j = 0 ; j < n ; j++ )
        {

            coeffMatrix[i][j] = coeffMatrixPhase1[i][j] ;

        }
    }
    //using coeffMatrix and phase1 basis vector in phase2
    displayCoeffiecientMatrix();
    cout<<endl;
    displayPhase1Basis();//using basis vector of phase 1

}

//make tableu legitimiate
    for(int i = 0; i < basisVectorPhase1.size(); i++)
    {
        //cout <<"x"<<basisVectorPhase1[i].first<<" = "<<basisVectorPhase1[i].second<<endl;
        int colofBasis = basisVectorPhase1[i].first - 1;
        if(coeffMatrix[i+1][colofBasis] == 1)
        {
            //make tableau a proper tableau
            int multplier = coeffMatrix[0][colofBasis];
            for(int k = 0; k < n; k++)//each column k of rowj effect row 0
            {
                coeffMatrix[0][k] = coeffMatrix[0][k] - ((multplier)*coeffMatrix[i+1][k]);
            }
        }
        else{
        cout<<"makin table lgitmate in phase 2 element not 1"<<endl;
        }
    }

    cout<<"legitimate tableau"<<endl;
    displayCoeffiecientMatrix();
//#if 0
    //now solve  within loop iterations
    while(1)
    {
        int entering_index =  searchMinz();
        if((entering_index == -1)||(coeffMatrix[0][entering_index] > 0))
        {
            cout<< "stopping simplex iterations" << endl;
            break;//optimal or unbounded solution
        }
        int leaving_index = calMinRatio_LeavingBasisVector(entering_index);

        if(leaving_index == -1)
        {
            cout<< "stopping simplex iterations" << endl;
            break;
        }

        double pivot_element = coeffMatrix[leaving_index][entering_index];

        //make pivot element unity along with row operation by dividing whole row by pivot element
        if(pivot_element != 0){
            for(int j = 0; j < n; j++)
            {
                coeffMatrix[leaving_index][j] /= pivot_element;
            }
        }
        else{
            cout<< "pivot element is 0"<< endl;
        }

        pivot_element = coeffMatrix[leaving_index][entering_index];//resetting the pivot element to the value after the roe operation ie 1
        cout<<"pivot = "<<pivot_element<<endl;
        cout<<"leaving row = "<<leaving_index<<" "<<"entering column = "<<entering_index<<endl;
        //make all elements of pivot column 0 expect pivot ele, Ri = Ri - Rp*Ri
        for( int i = 0 ; i < m ; i++ )
        {
            double multiply_factor = coeffMatrix[i][entering_index];
            if(leaving_index != i)
            {
                for(int j = 0; j < n; j++)
                {
                    coeffMatrix[i][j] = coeffMatrix[i][j] - ((multiply_factor) * (coeffMatrix[leaving_index][j]));
                    if(fabs(coeffMatrix[i][j]) < 0.000001)
                    {
                        cout<< "setting 0  "<<coeffMatrix[i][j]<<endl;
                        coeffMatrix[i][j] = 0;

                    }

                }

            }
        }

        updateSolutionvectorPhase2(leaving_index, entering_index);
        cout <<"simplex tableau at current iteration : "<<endl;

        displayCoeffiecientMatrix();
        cout<<endl;
        displayPhase1Basis();//note: using phase1 basis vector
        cout<<endl<<"********************************************************************************"<<endl;
    }

//#endif
}

/*
desp :  finds if the colIndex is on amongst the artificial variables
author : Ajay Singh
*/
bool artColinBasis(int colIndex)//col index starts from 1
{
    for(int i = 0; i < basisVectorPhase1.size(); i++)
    {
        if(basisVectorPhase1[i].first == colIndex) //belongs to artificial vars col
        {
            return true;
        }
    }
    return false;
}

