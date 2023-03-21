#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <cmath>

using namespace std;

const double eps = 1e-12;


struct Bok {
    int iloscPC; 
    double* wspKsi;
    double* wspEta;
    double** wartFK; 
   
};

struct Elem_Uni {
    double** Xi_tab;
    double** N_tab;
    int L_puntk;
    Bok* boki;
};


struct Pc_Mat {
    double** NX_tab;
    double** NY_tab;
    double* PcJacob;
    int L_puntk;
    Bok* boki;
};

struct Node {
    double x;
    double y;
    double t;
    int BC;
};

struct Element {
    int ID[4]; 
};

struct Grid {
    int nE;
    int nN;
    Node* ND;
    Element* NE;

};

struct H_grid {
    double** H_mat;
    double** C_mat;
    double* P_vec;
    int L_punkt;
};

struct GlobalData {
    int alpha;
    int t_otoczenia; //tot
    int cp; //specificHeat
    int density; //density
    int t; //simulation time
    int delta_t; //simulationStepTime
    int conductivity;
    int initialTemp;

};


void getData(GlobalData* gd, Grid* grid);
string* splitstr(string str, int n, char deli);
int str_to_int(string str);
string delSpace(string str);
int countIn(string str, char deli);

//lab2  
double F_x(double x);
double F_nj2D(double e, double n);
double calka(int n, double* PC, double* W, int przestrzen);

//lab3
double calkaJakobian(double* PC, double* W, int point, double a, double b);
double F_x2(double x);


//lab4
double N_xi(int LP, double n);
double N_et(int LP, double x);
void parseToElemUni(Elem_Uni* uni, int L_punkt, double pktTabXi[], double pktTabEt[]); //tworzenie elementu uniwersalnego 
void ShowUni(Elem_Uni uni, double pktTabXi[], double pktTabEt[]);

//lab5
void parseToPcMat(Pc_Mat* _pcmat, int L_punkt, Elem_Uni _uni, double X_tab[], double Ytab[], Grid _grid); //tworzy macierz dN/dx i dN/dY//_pcmat - tablica  dN/dx(dy) //L_punkt - PC //_uni - dN/dksi(deta) //punkty X //punkty Y
void parseToH(H_grid* _hgr, Pc_Mat _pcm, double w1, double w2, double Xtab[], double Ytab[]); //czysta macierz lokalna H bez BC
void ShowNxNy(Pc_Mat pc);

//lab6 AGREGACJA
void Agregacja(H_grid* Ha, H_grid hgr, Grid grid, Pc_Mat pcMat);

//lab7
double* NMat(double ksi, double eta);

//lab8 TEMPERATURA ZE WZORU Ht + P = 0
void Temp(H_grid Ha, Grid grid);

//lab9 TEMPERATURA Z KROKIEM CZASOWYM
void stepTempTime(Grid grid, H_grid hgr);

//funkcje do ostatniego równania
double minFunction(double* T, int size);
double maxFunction(double* T, int size);
double* GaussElimination(double** HG, double* P, int nN);



//zmienne globalne
GlobalData global_data;

int main() {
    ofstream zapis1("macierzHBC.txt");
    ofstream zapis2("macierzC.txt");
    ofstream zapis3("wektorP.txt");
   

    //----------lab1------------------------------------------
    Grid grid;
    string lin = "";

    getData(&global_data, &grid);
    

    //-------------lab2-----------------------------------------
    cout.precision(11);
    int n1 = 1;
    int n2 = 2;

    //punkty calkowania dla 2pkt i 3 pkt
    double PC2[2] = { -1.0 / sqrt(3.0),1.0 / sqrt(3.0) };
    double PC3[3] = { -sqrt(3.0 / 5.0),0.0,sqrt(3.0 / 5.0) };

    //wagi dla 2pkt i 3pkt
    double W1[2] = { 1.0,1.0 };
    double W2[3] = { double(5.0 / 9.0),double(8.0 / 9.0),double(5.0 / 9.0) };


    double C1 = calka(n1, PC2, W1, 1);
    double C1_2 = calka(n2, PC3, W2, 1);

    double C2 = calka(n1, PC2, W1, 2);
    double C2_2 = calka(n2, PC3, W2, 2);


    //--------------lab3------------------------------------------
    double value = calkaJakobian(PC3, W2, 3, 3, 8);
   
    //------------lab4-----------------PROJEKT-------------------
    Elem_Uni _uni;
    double XiWart[4] = { -0.57735,0.57735,0.57735,-0.57735 }; //wspolrzedne punktow calkowania
    double EtWart[4] = { -0.57735,-0.57735,0.57735,0.57735 };
    int points = 4; //2d //4pc //liczba punktow calkowania

    parseToElemUni(&_uni, points, XiWart, EtWart); //tworzenie el uniwersalnego
    ShowUni(_uni, XiWart, EtWart); //wypisywanie wartosci elementu


    //-----------lab5-------------------------------------------

    Pc_Mat pc;
    double _Xtab[4]; //={ 0.0,0.025,0.025,0.0 };  
    double _Ytab[4]; //= { 0.0,0.0,0.025,0.025 };

    for (int i = 0; i < 4; i++) {
        _Xtab[i] = grid.ND[grid.NE[0].ID[i]].x;
        _Ytab[i] = grid.ND[grid.NE[0].ID[i]].y;
    }

    H_grid hgr;

    parseToPcMat(&pc, points, _uni, _Xtab, _Ytab, grid); // macierz: dN/dx //dN/dy

    
    ShowNxNy(pc); //wypisanie
    cout << endl << endl;

    parseToH(&hgr, pc, 1.0, 1.0, XiWart, EtWart); //lokalna macierz H

    //macierz H lokalna
    cout << "Macierz H lokalna:" << endl;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            cout << hgr.H_mat[i][j] << "\t";
        }
        cout << endl;
    }

    H_grid H_agr; 

    //--------------lab6

    Agregacja(&H_agr, hgr, grid, pc);


    cout << endl << endl;
    cout << "Macierz H+HBC:" << endl;
    zapis1 << "Macierz H+HBC:" << endl;
    for (int i = 0; i < grid.nN; i++) {
        for (int j = 0; j < grid.nN; j++) {
            cout << std::setprecision(4) << H_agr.H_mat[i][j] << "\t";
            zapis1 << std::setprecision(4) << H_agr.H_mat[i][j] << "\t";
        }
        cout << endl;
        zapis1 << endl;
    }


    cout << endl << endl;
    cout << "Wektor P:" << endl;
    zapis3 << "Wektor P:" << endl;
    for (int i = 0; i < grid.nN; i++) {
        cout << i + 1 << ": " << setprecision(8) << H_agr.P_vec[i] << endl;
        zapis3 << i + 1 << ": " << setprecision(8) << H_agr.P_vec[i] << endl;
    }

    //---------------lab7

    cout << endl << endl;

    //TEMPERATURA Z Ht+P=0
    //Temp(H_agr, grid);


    cout << endl << endl;

    cout << "Macierz C:" << endl << endl;
    zapis2 << "Macierz C:" << endl << endl;
    for (int i = 0; i < grid.nN; i++) {
        for (int j = 0; j < grid.nN; j++) {
            cout << H_agr.C_mat[i][j] << "\t";
            zapis2 << H_agr.C_mat[i][j] << "\t";
        }
        cout << endl;
        zapis2 << endl;
    }


    cout << endl << endl << endl;

    //TEMPERATURA wyliczona z ostatniego rownania
    stepTempTime(grid, H_agr);



	zapis1.close();
	zapis2.close();
	zapis3.close();
}

string* splitstr(string str, int n, char deli) //linia, ile rzeczy do podzielenia, znak podzialu
{
    string* Tab = new string[n];
    int counter = 0;
    for (int i = 0; i < str.size(); i++) {

        if (str[i] == deli) {
            counter++;
            continue;
        }

        Tab[counter] = Tab[counter] + str[i];
    }
    return Tab;

}

int str_to_int(string str)
{
    stringstream ss;
    int num;

    ss << str;
    ss >> num;

    return num;

}


string delSpace(string str)
{
    string s = "";
    for (int i = 0; i < str.size(); i++) {
        if (str[i] == ' ') continue;
        else s = s + str[i];
    }
    return s;
}

int countIn(string str, char deli) {
    int counter = 0;
    for (int i = 0; i < str.size(); i++) {
        if (str[i] == deli) counter++;
    }
    return counter;
}



void getData(GlobalData* gd, Grid* grid) {
    int i = 1;
    //ifstream odczyt("siatka_dodatkowa.txt");
    ifstream odczyt("Test1_4_4.txt");
    //ifstream odczyt("Test2_4_4_MixGrid.txt");
    //ifstream odczyt("Test3_31_31_kwadrat.txt"); 

   

 
    gd->t = 1.0;
    string linia;

    Node* NodeTab = new Node[0];
    Element* ElementTab = new Element[0];


    while (getline(odczyt, linia)) {

        if (linia[0] == '*') {

            if (linia[1] == 'N') {
                i = 3;
                NodeTab = new Node[grid->nN];	
                continue;
            }
            else if (linia[1] == 'E') {
                grid->ND = NodeTab;
                ElementTab = new Element[grid->nE];

                i = 4;
                continue;
            }
            else if (linia[1] == 'B') {
                grid->NE = ElementTab;
                i = 5;
                continue;
            }
        }
        if (i == 1) {
            string* tablica = new string[2];
            tablica = splitstr(linia, 2, ' ');
            if (tablica[0] == "SimulationTime") {
                gd->t = str_to_int(tablica[1]);

            }
            else if (tablica[0] == "SimulationStepTime") {
                gd->delta_t = str_to_int(tablica[1]);
            }
            else if (tablica[0] == "Conductivity") {
                gd->conductivity = str_to_int(tablica[1]);
            }
            else if (tablica[0] == "Alfa") {
                gd->alpha = str_to_int(tablica[1]);
            }
            else if (tablica[0] == "Tot") {
                gd->t_otoczenia = str_to_int(tablica[1]);
            }
            else if (tablica[0] == "InitialTemp") {
                gd->initialTemp = str_to_int(tablica[1]);
            }
            else if (tablica[0] == "Density") {
                gd->density = str_to_int(tablica[1]);
            }
            else if (tablica[0] == "SpecificHeat") {
                gd->cp = str_to_int(tablica[1]);


                i = 2;
            }


            //string *parts = linia.splitstr(" ", 2);
            //dzielenie na 2 czesci, do tablicy 2 elementowej
            //zamiana drugiego elementu na wart liczbowa
            //(zamiana na int, ew. zamiana wszystkich double na int)
            //2. porownanie pierwszego elementu z tablicy z wartosciami ze struktury


        }
        else if (i == 2) {
            string* tablica2 = new string[3];
            tablica2 = splitstr(linia, 3, ' ');
            if (tablica2[0] == "Nodes") {
                grid->nN = str_to_int(tablica2[2]);
            }
            else if (tablica2[0] == "Elements") {
                grid->nE = str_to_int(tablica2[2]);
            }
        }
        else if (i == 3) {
            string linia2 = delSpace(linia);
            string* tabNodes = new string[3];
            tabNodes = splitstr(linia2, 3, ',');

            //zamiana na double
            double x = stod(tabNodes[1]);
            double y = stod(tabNodes[2]);


            int licznik = str_to_int(tabNodes[0]);
            licznik -= 1;
            NodeTab[licznik].x = x;
            NodeTab[licznik].y = y;
            NodeTab[licznik].t = gd->initialTemp;
            NodeTab[licznik].BC = 0; //ustawienie 0 dla wszystkich 
           
            //podzial stringa na 3 elementy po [, spacja]
            //funk. tworzenie node, przypisywanie wspolrz.

        }
        else if (i == 4) {

            string linia3 = delSpace(linia);
            string* tabElements = new string[5];
            tabElements = splitstr(linia3, 5, ',');

            for (int k = 0; k < 4; k++) {
                ElementTab[str_to_int(tabElements[0]) - 1].ID[k] = str_to_int(tabElements[k + 1]);
            }



            //podzial stringa na 5 elementow podobnie jak w poprzednim
            //funkt. tworzenie elementow, przypisywanie wczesniej stworzonych node
        }
        else if (i == 5) {
            string linia4 = delSpace(linia);
            int liczbaPrzecinkow = countIn(linia4, ',');

            string* tabBC = new string[liczbaPrzecinkow + 1];

            tabBC = splitstr(linia4, liczbaPrzecinkow + 1, ',');

            for (int j = 0; j <= liczbaPrzecinkow; j++) {
                int index = str_to_int(tabBC[j]) - 1;
                grid->ND[index].BC = 1; //usatwienie 1 dla wybranych
            }
           

            //a w 4 na n elementow, bazujac na ilosci przecinkow i 
            //sprawdzajac czt nie nie skonczyla linia przypadkiem
            //funk.modyfikacja konkretnych node (zamiana flagi na brzegowa)
            break;
        }

        //if zczytujacy pierwsze dwa chary z linii, jezeli one beda rowne *N to zamieniamy i na 2
        //a jezeli beda rowne *E to zamieniamy i na 3, a jesli *B to na 4
    }

    odczyt.close();





}

//----------lab2-----------------------------

double F_x(double x)
{
    return ((2 * x * x) + (3 * x) - 8);
}

double F_nj2D(double e, double n)
{
    return -5 * e * e * n + 5 * e * n * n + 10;
}

double calka(int n, double* PC, double* W, int przestrzen)
{
    double suma = 0;
    if (przestrzen == 1)
    {
        for (int i = 0; i <= n; i++)
        {
            suma += F_x(PC[i]) * W[i];
        }
    }
    else if (przestrzen == 2)
    {
        for (int i = 0; i <= n; i++)
        {
            for (int j = 0; j <= n; j++)
            {
                suma += F_nj2D(PC[i], PC[j]) * W[i] * W[j];
            }
        }
    }
    else
    {
        cout << "Przestrzen bledna" << endl;

    }
    return suma;
}

//-----------lab3-------------------------------
double F_x2(double x)
{
    return ((x * x) - 3 * x + 6);
}

//FUNKCJE KSZTALTU
//funkcja liczaca pochodna w od wartosci (n) w zaleznosci od funkcji N(LP)
double N_xi(int LP, double n) 
{
    switch (LP) {
    case 1:
        return (double)(-0.25 * (1.0 - n));
        break;
    case 2:
        return (double)(0.25 * (1.0 - n));
        break;

    case 3:
        return (double)(0.25 * (1.0 + n));
        break;

    case 4:
        return (double)(-0.25 * (1.0 + n));
        break;
    default:
        return 0.0;
        break;
    }
}

//funkcja liczaca pochodna w od wartosci (n) w zaleznosci od funkcji N(LP)
double N_et(int LP, double x) 
{
    switch (LP) {
    case 1:
        return (double)(-0.25 * (1.0 - x));
        break;
    case 2:
        return (double)(-0.25 * (1.0 + x));
        break;

    case 3:
        return (double)(0.25 * (1.0 + x));
        break;

    case 4:
        return (double)(0.25 * (1.0 - x));
        break;

    default:
        return 0.0;
        break;

    }
}


void parseToElemUni(Elem_Uni* uni, int L_punkt, double pktTabXi[], double pktTabEt[])
{
    //tworzenie pustych wartosci
    double** XiTab = new double* [L_punkt];
    for (int i = 0; i < 4; i++) {
        XiTab[i] = new double[4];
    }

    //tworzenie pustych wartosci
    double** EtTab = new double* [L_punkt];
    for (int i = 0; i < 4; i++) {
        EtTab[i] = new double[4];
    }

    //libczby punktow calkowania - dla psi
    for (int i = 0; i < L_punkt; i++) {
        for (int j = 0; j < 4; j++) {
            XiTab[i][j] = N_xi(j + 1, pktTabEt[i]);
        }
    }

    //liczby punktow calkowania - dla eta
    for (int i = 0; i < L_punkt; i++) {
        for (int j = 0; j < 4; j++) {
            EtTab[i][j] = N_et(j + 1, pktTabXi[i]);
        }
    }

    Bok* boki = new Bok[4];

    for (int i = 0; i < 4; i++) {
        boki[i].iloscPC = 2; //pc na boku //kazdy ma 2
    }

    //wspolrzedne dla PC na boku //po 2 PC na bok 
    for (int i = 0; i < 4; i++) {
        boki[i].wspKsi = new double[2];
        boki[i].wspEta = new double[2];
    }

    //Tworzenie bokow

    //bok1
    boki[0].wspKsi[0] = -0.5773;
    boki[0].wspKsi[1] = 0.5773;
    boki[0].wspEta[0] = -1.0;
    boki[0].wspEta[1] = -1.0;

    //bok2
    boki[1].wspKsi[0] = 1.0;
    boki[1].wspKsi[1] = 1.0;
    boki[1].wspEta[0] = -0.5773;
    boki[1].wspEta[1] = 0.5773;

    //bok3
    boki[2].wspKsi[0] = 0.5773;
    boki[2].wspKsi[1] = -0.5773;
    boki[2].wspEta[0] = 1.0;
    boki[2].wspEta[1] = 1.0;

    //bok4
    boki[3].wspKsi[0] = -1.0;
    boki[3].wspKsi[1] = -1.0;
    boki[3].wspEta[0] = 0.5773;
    boki[3].wspEta[1] = -0.5773;


    //wartosci N dla pc od ksi i eta //ta macierz kolo tego rysunku z pc
    for (int i = 0; i < 4; i++) {
        boki[i].wartFK = new double* [2];
        for (int j = 0; j < 2; j++) {
            boki[i].wartFK[j] = new double[4];
            for (int z = 0; z < 4; z++) {
                boki[i].wartFK[j][z] = NMat(boki[i].wspKsi[j], boki[i].wspEta[j])[z];
            }


            cout << boki[i].wspKsi[j] << "\t" << boki[i].wspEta[j] << "\t";
            for (int z = 0; z < 4; z++) {
                cout << boki[i].wartFK[j][z] << "\t";
            }
            cout << endl;

        }



    }

    uni->L_puntk = L_punkt;
    uni->N_tab = EtTab;
    uni->Xi_tab = XiTab;
    uni->boki = boki;
}

//------lab5
void parseToPcMat(Pc_Mat* _pcmat, int L_punkt, Elem_Uni _uni, double Xtab[], double Ytab[], Grid _grid) {

    _pcmat->boki = _uni.boki;

    //tworzenie pustych wartosci
    double** NxTab = new double* [L_punkt];
    for (int i = 0; i < 4; i++) {
        NxTab[i] = new double[4];
    }

    //tworzenie pustych wartosci
    double** NyTab = new double* [L_punkt];
    for (int i = 0; i < 4; i++) {
        NyTab[i] = new double[4];
    }

    //rozmiar tablicy Jakobian dla pc
    double* PcJacobTab = new double[L_punkt];

    for (int pc = 0; pc < L_punkt; pc++) { // dla kazcego punktu PC



        double DyDet, DyDxi, DxDet, DxDxi;

        DyDet = 0.0;//
        for (int i = 0; i < 4; i++) {

            DyDet += _uni.N_tab[pc][i] * Ytab[i];
        }//policzone


        DyDxi = 0.0;//
        for (int i = 0; i < 4; i++) {
            DyDxi += _uni.Xi_tab[pc][i] * Ytab[i];

        }//policzone


        DxDet = 0.0;//
        for (int i = 0; i < 4; i++) {

            DxDet += _uni.N_tab[pc][i] * Xtab[i];
        }//policzone


        DxDxi = 0.0;//
        for (int i = 0; i < 4; i++) {
            DxDxi += _uni.Xi_tab[pc][i] * Xtab[i];
        }//policzone


        //gotowa macierz ta dx po det dla konkretnego punktu

        double det = (DxDxi * DyDet) - (DxDet * DyDxi);
        PcJacobTab[pc] = det;                               // wpisanie Jacobianu dla konkretnego Punktu calkowania (PC)
        det = 1.0 / det;
        DxDxi *= det;
        DxDet *= det;
        DyDxi *= det;
        DyDet *= det;


        for (int i = 0; i < 4; i++) {
            NxTab[pc][i] = DxDxi * _uni.Xi_tab[pc][i] + DyDxi * _uni.N_tab[pc][i];
            NyTab[pc][i] = DxDet * _uni.Xi_tab[pc][i] + DyDet * _uni.N_tab[pc][i];
        }

    }

    //wkladamy wartosci do struktury pcmat

    _pcmat->L_puntk = L_punkt;
    _pcmat->NX_tab = NxTab;
    _pcmat->NY_tab = NyTab;

    _pcmat->PcJacob = PcJacobTab;

}

void parseToH(H_grid* _hgr, Pc_Mat _pcm, double w1, double w2, double Xtab[], double Ytab[])
{
    //tworzenie pustych macierzy H dla poszczegolnych punktow calkowania
    // 
    //macierze H dla kazdego PC
    //macierz C dla kazdego PC
    double*** Hx = new double** [_pcm.L_puntk]; 
    double*** Cx = new double** [_pcm.L_puntk];
    
    for (int i = 0; i < 4; i++) {
        Hx[i] = new double* [4];
        Cx[i] = new double* [4];
        for (int j = 0; j < 4; j++) {
            Hx[i][j] = new double[4];
            Cx[i][j] = new double[4];
        }
    }

    double*** Hy = new double** [_pcm.L_puntk];
    double*** Cy = new double** [_pcm.L_puntk];
    for (int i = 0; i < 4; i++) {
        Hy[i] = new double* [4];
        Cy[i] = new double* [4];
        for (int j = 0; j < 4; j++) {
            Hy[i][j] = new double[4];
            Cy[i][j] = new double[4];
        }
    }

    double*** H_lok = new double** [_pcm.L_puntk];
    double*** C_lok = new double** [_pcm.L_puntk];
    for (int i = 0; i < 4; i++) {
        H_lok[i] = new double* [4]; //macierz lokalna H dla wszystkich PC
        C_lok[i] = new double* [4]; //macierz lokalna C dla wszystkich PC
        for (int j = 0; j < 4; j++) {
            H_lok[i][j] = new double[4];
            C_lok[i][j] = new double[4];
        }
    }

    for (int p = 0; p < _pcm.L_puntk; p++) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                C_lok[p][i][j] = 0.0;
            }
        }
    }

    double** HH = new double* [4]; //macierz lokalna H dla jednego PC
    double** CC = new double* [4]; //macierz lokalna C dla jednego PC
    for (int i = 0; i < 4; i++) {
        HH[i] = new double[4];
        CC[i] = new double[4];
    }

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            HH[i][j] = 0.0;
            CC[i][j] = 0.0;
        }
    }

    for (int pc = 0; pc < _pcm.L_puntk; pc++) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {

                // tworzenie macierzy co wchodza do rownania H dla kazdego punktu calkowania
                Hx[pc][i][j] = _pcm.NX_tab[pc][i] * _pcm.NX_tab[pc][j]; //pierwsza
                Hy[pc][i][j] = _pcm.NY_tab[pc][i] * _pcm.NY_tab[pc][j]; //druga


                Cx[pc][i][j] = _pcm.NX_tab[pc][i] * _pcm.NX_tab[pc][j]; //pierwsza
                Cy[pc][i][j] = _pcm.NY_tab[pc][i] * _pcm.NY_tab[pc][j]; //druga


            }
        }

        // ten wyznacznik H dla kazdego punktu calkowania

//macierz H przemnozona przez jacobian
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
              
                //dla wszystkich pc
                H_lok[pc][i][j] = global_data.conductivity * _pcm.PcJacob[pc] * (Hx[pc][i][j] + Hy[pc][i][j]); 
                C_lok[pc][i][j] = global_data.cp * global_data.density * _pcm.PcJacob[pc] * NMat(Xtab[pc], Ytab[pc])[i] * NMat(Xtab[pc], Ytab[pc])[j]; //(NMat(Xtab[pc], Ytab[pc])[i] * NMat(Xtab[pc], Ytab[pc])[j]); 
                
               
            }
        }
    }

    //Ostatnie obliczenia - wynik to macierz H
    for (int pc = 0; pc < _pcm.L_puntk; pc++) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                                                        //cout << H_lok[pc][i][j] << endl;
                HH[i][j] += H_lok[pc][i][j] * w1 * w2; //sumujemy dla wszystkich pc
                CC[i][j] += C_lok[pc][i][j] * w1 * w2; //sumujemy dla wszystkich pc
            }
        }
    }

    _hgr->L_punkt = _pcm.L_puntk;
    _hgr->H_mat = HH;
    _hgr->C_mat = CC;

}

void ShowNxNy(Pc_Mat pc)
{
    cout << endl << endl;
    cout << "Nx" << endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            cout << pc.NX_tab[i][j] << "\t";
        }
        cout << endl;
    }

    cout << endl;
    cout << "Ny" << endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            cout << pc.NY_tab[i][j] << "\t";
        }
        cout << endl;
    }

    //jakobian dla PC
    cout << endl;
    for (int i = 0; i < pc.L_puntk; i++) {
        cout << "Jacobian" << i + 1 << ": " << pc.PcJacob[i] << endl;
    }

}


//wyswietlenie macierzy ksi i eta
void ShowUni(Elem_Uni uni, double pktTabXi[], double pktTabEt[])
{
    cout << endl << endl;
    cout << ".\teta\tdN1/dksi\tdN2/dksi\tdN3/dksi\tdN4/dksi" << endl;
    for (int i = 0; i < uni.L_puntk; i++) {
        cout << "PC" << i + 1 << "\t";
        cout << pktTabEt[i] << "\t";
        for (int j = 0; j < 4; j++) {
            cout << uni.Xi_tab[i][j] << "\t";
        }
        cout << endl;
    }

    cout << endl << endl;
    cout << ".\tksi\tdN1/deta\tdN2/deta\tdN3/deta\tdN4/deta" << endl;
    for (int i = 0; i < uni.L_puntk; i++) {
        cout << "PC" << i + 1 << "\t";
        cout << pktTabXi[i] << "\t";
        for (int j = 0; j < 4; j++) {
            cout << uni.N_tab[i][j] << "\t";
        }
        cout << endl;
    }

}

double calkaJakobian(double* PC, double* W, int point, double a, double b) {
    double value = 0;
    double detJ = (b - a) / 2;
    for (int i = 0; i < point; i++)
    {
        PC[i] = ((1 - PC[i]) / 2) * a + ((PC[i] + 1) / 2) * b;
        value += (F_x2(PC[i]) * W[i]);
    }

    value = value * detJ;
    return value;
}


//TODO ogarnac tak zeby dzialalo
void Agregacja(H_grid* Ha, H_grid hgr, Grid grid, Pc_Mat pcMat) {
    Ha->L_punkt = hgr.L_punkt;


    cout << endl << endl;
    for (int i = 0; i < grid.nN; i++) {
        cout << i << "\t" << grid.ND[i].BC << endl;
    }

    cout << endl << endl;

    double* Pe_V = new double[grid.nN];
    double** Ha_Mat = new double* [grid.nN];
    double** Ca_Mat = new double* [grid.nN];
    for (int i = 0; i < grid.nN; i++) {
        Ha_Mat[i] = new double[grid.nN];
        Ca_Mat[i] = new double[grid.nN];
    }
    for (int i = 0; i < grid.nN; i++) {
        Pe_V[i] = 0.0;
        for (int j = 0; j < grid.nN; j++) {
            Ha_Mat[i][j] = 0.0;
            Ca_Mat[i][j] = 0.0;


        }
    }

    //pkt w elemncie
    for (int pkt = 0; pkt < grid.nE; pkt++) {
        for (int i = 0; i < 4; i++) {
            cout << grid.NE[pkt].ID[i] << "\t";
        }


        cout << endl << endl;

        double** HbcTab1 = new double* [hgr.L_punkt];
        double* PTab1 = new double[4];

        for (int i = 0; i < 4; i++) {
            HbcTab1[i] = new double[4];
        }

        double** Hbc = new double* [hgr.L_punkt];
        double* PTabV = new double[4];

        for (int i = 0; i < 4; i++) {
            Hbc[i] = new double[4];
        }

        //tworzenie pustych wartosci
        double** HbcTab2 = new double* [hgr.L_punkt];
        double* PTab2 = new double[4];
        for (int i = 0; i < 4; i++) {
            HbcTab2[i] = new double[4];
        }

        for (int i = 0; i < 4; i++) {
            PTabV[i] = 0.0;
            for (int j = 0; j < 4; j++) {
                Hbc[i][j] = 0.0;
            }
        }

        for (int bokI = 0; bokI < 4; bokI++) {


            for (int i = 0; i < 4; i++) {

                //PTab1
                PTab1[i] = pcMat.boki[bokI].wartFK[0][i] * global_data.t_otoczenia * global_data.alpha;  //1Pc
                //PTab2
                PTab2[i] = pcMat.boki[bokI].wartFK[1][i] * global_data.t_otoczenia * global_data.alpha; //2Pc

                for (int j = 0; j < 4; j++) {
                    //HBC1
                    HbcTab1[i][j] = pcMat.boki[bokI].wartFK[0][i] * pcMat.boki[bokI].wartFK[0][j] * global_data.alpha;
                    //HBC2
                    HbcTab2[i][j] = pcMat.boki[bokI].wartFK[1][i] * pcMat.boki[bokI].wartFK[1][j] * global_data.alpha;

                }
            }
           

            int nodeNum1 = grid.NE[pkt].ID[bokI] - 1;
            int nodeNum2;

            //szukanie id node dla boku
            if (bokI == 3) {
                nodeNum2 = grid.NE[pkt].ID[0] - 1;
            }
            else {
                nodeNum2 = grid.NE[pkt].ID[bokI + 1] - 1;
            }

            double dlBoku = sqrt(pow(grid.ND[nodeNum1].x - grid.ND[nodeNum2].x, 2) + pow(grid.ND[nodeNum1].y - grid.ND[nodeNum2].y, 2));

            double detJboku = dlBoku / 2.0;

           /* cout << endl << endl << "COND: " << global_data.conductivity << " dETJ: " << detJboku << "\t DlBoku: " << dlBoku << endl;

            cout << endl << endl;
            cout << "nodeNum1: " << nodeNum1 << " BC: " << grid.ND[nodeNum1].BC << "\tnodeNum2: " << nodeNum2 << " BC:" << grid.ND[nodeNum2].BC << endl;*/

            for (int i = 0; i < 4; i++) {

                //liczenie wektora P
                PTabV[i] = PTabV[i] + ((PTab1[i] + PTab2[i]) * detJboku * (grid.ND[nodeNum1].BC * grid.ND[nodeNum2].BC)); //sprytny plan

                for (int j = 0; j < 4; j++) {
                    Hbc[i][j] = Hbc[i][j] + ((HbcTab1[i][j] + HbcTab2[i][j]) * detJboku * (grid.ND[nodeNum1].BC * grid.ND[nodeNum2].BC));
                }
            }



        }

        for (int i = 0; i < 4; i++) {
            Pe_V[grid.NE[pkt].ID[i] - 1] += PTabV[i]; //pev - glob //ptab -lok
            for (int j = 0; j < 4; j++) {
                // Hbc[i][j] = Hbc[i][j] + ((HbcTab1[i][j] + HbcTab2[i][j]) * detJboku * (grid.ND[nodeNum1].BC * grid.ND[nodeNum2].BC));
                Ha_Mat[grid.NE[pkt].ID[i] - 1][grid.NE[pkt].ID[j] - 1] += Hbc[i][j];

            }
        }

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                cout << setprecision(2) << Hbc[i][j] << "\t";
            }
            cout << endl;
        }
        cout << endl << endl;


        for (int i = 0; i < 4; i++) {

            cout << setprecision(6) << PTabV[i] << "\t";
        }
        cout << endl << endl;



    }

    //AgregatioMes

    for (int pkt = 0; pkt < grid.nE; pkt++) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {

                Ha_Mat[grid.NE[pkt].ID[i] - 1][grid.NE[pkt].ID[j] - 1] += hgr.H_mat[i][j]; 
                Ca_Mat[grid.NE[pkt].ID[i] - 1][grid.NE[pkt].ID[j] - 1] += hgr.C_mat[i][j];
            }
        }
    }


    Ha->P_vec = Pe_V;
    Ha->H_mat = Ha_Mat;
    Ha->C_mat = Ca_Mat;


}

double* NMat(double ksi, double eta)
{
    double N1 = 0.25 * (1.0 - ksi) * (1.0 - eta);
    double N2 = 0.25 * (1.0 + ksi) * (1.0 - eta);
    double N3 = 0.25 * (1.0 + ksi) * (1.0 + eta);
    double N4 = 0.25 * (1.0 - ksi) * (1.0 + eta);


    double mat[4] = { N1,N2,N3,N4 };
    return mat;
}

bool ludist(int n, double** A)
{
    int i, j, k;

    for (k = 0; k < n - 1; k++)
    {
        if (fabs(A[k][k]) < eps) return false;

        for (i = k + 1; i < n; i++)
            A[i][k] /= A[k][k];

        for (i = k + 1; i < n; i++)
            for (j = k + 1; j < n; j++)
                A[i][j] -= A[i][k] * A[k][j];
    }

    return true;
}

// Funkcja wyznacza wektor X na podstawie A i X [ i ] 
//---------------------------------------------------
bool lusolve(int k, int n, double** A, double** X)
{
    int    i, j;
    double s;

    for (i = 1; i < n; i++)
    {
        s = 0;

        for (j = 0; j < i; j++) s += A[i][j] * X[j][k];

        X[i][k] -= s;
    }

    if (fabs(A[n - 1][n - 1]) < eps) return false;

    X[n - 1][k] /= A[n - 1][n - 1];

    for (i = n - 2; i >= 0; i--)
    {
        s = 0;

        for (j = i + 1; j < n; j++) s += A[i][j] * X[j][k];

        if (fabs(A[i][i]) < eps) return false;

        X[i][k] = (X[i][k] - s) / A[i][i];
    }

    return true;
}


void Temp(H_grid Ha, Grid grid)
{
    double** A, ** X;
    int i, j;
    bool ok;


    A = new double* [grid.nN];
    X = new double* [grid.nN];

    for (i = 0; i < grid.nN; i++)
    {
        A[i] = new double[grid.nN];
        X[i] = new double[grid.nN];
    }

    //zapisa danych dla macierzy A
    for (i = 0; i < grid.nN; i++)
        for (j = 0; j < grid.nN; j++)  A[i][j] = Ha.H_mat[i][j];


    if (ludist(grid.nN, A))
    {

        // Tworzymy macierz jednostkową w X    

        for (i = 0; i < grid.nN; i++)
        {
            for (j = 0; j < grid.nN; j++) X[i][j] = 0;
            X[i][i] = 1;
        }

        // Wyznaczamy kolejne kolumny macierzy odwrotnej w X

        ok = true;
        for (i = 0; i < grid.nN; i++)
            if (!lusolve(i, grid.nN, A, X))
            {
                ok = false;
                break;
            }
    }
    else ok = false;

    // Jeśli obliczenia się powiodły, wyświetlamy X

    cout << endl << endl;

  
    double* Tep_v = new double[grid.nN];

    for (int i = 0; i < grid.nN; i++) {
        Tep_v[i] = 0;
    }


    for (int i = 0; i < grid.nN; i++) {
        for (j = 0; j < grid.nN; j++) {
            Tep_v[i] += X[i][j] * (-1.0 * Ha.P_vec[i]);  //rozwiazanie równ tep
        }
    }

    cout << endl << endl;
    cout << "Temp z rownania Ht+P=0" << endl;
    for (int i = 0; i < grid.nN; i++) {

        cout << Tep_v[i]<< endl;


    }


}

double minFunction(double* T, int size) {
    double min = T[0];
    for (int i = 0; i < size; i++)
    {
        if (min > T[i]) min = T[i];
    }
    return min;
}

double maxFunction(double* T, int size) {
    double max = T[0];
    for (int i = 0; i < size; i++)
    {
        if (max < T[i]) max = T[i];
    }
    return max;
}

//-------------Gauss
double* GaussElimination(double** HG, double* P, int nN)
{
    double** result = new double* [nN];
    double* temperatura = new double[nN];
   
    int i, j, k;
    double m, s;
    const double eps = 1e-12;

    for (int i = 0; i < nN; i++) result[i] = new double[nN + 1];

    for (int i = 0; i < nN; i++) {
        for (int j = 0; j < nN; j++)
            result[i][j] = HG[i][j];
    }

    for (int i = 0; i < nN; i++) result[i][nN] = P[i];
    
    for (int i = 0; i < nN; i++) temperatura[i] = 0.0;
    
    for (i = 0; i < nN - 1; i++) {
        for (j = i + 1; j < nN; j++) {
            if (fabs(result[i][i]) < eps) break;
                m = -result[j][i] / result[i][i];
            for (k = i + 1; k <= nN; k++)
                result[j][k] += m * result[i][k];
        }
    }
    for (i = nN - 1; i >= 0; i--) {
        s = result[i][nN];
        for (j = nN - 1; j >= i + 1; j--)
            s -= result[i][j] * temperatura[j];
        if (fabs(result[i][i]) < eps) break;
        temperatura[i] = s / result[i][i];
    }
    return temperatura;
}



//-----------------------WYLICZENIE TEMPERATURY
void stepTempTime(Grid grid, H_grid hgr) {
   // ofstream zapis4("temp.txt");

    double tempMin, tempMax;
    double** sumaHC = new double* [grid.nN]; //suma macierzy H i macierzy C
    double* sumaCP = new double(grid.nN); //suma C i wektora P
    double* t0 = new double(grid.nN);

    for (int i = 0; i < grid.nN; i++) sumaHC[i] = new double[grid.nN];
    for (int i = 0; i < grid.nN; i++) t0[i] = global_data.initialTemp;

    for (int i = 0; i < grid.nN; i++) {
        for (int j = 0; j < grid.nN; j++) {
            sumaHC[i][j] = (hgr.C_mat[i][j] / global_data.delta_t) + hgr.H_mat[i][j];
        }
    }

    for (int step = 0; step < global_data.t; step += global_data.delta_t) {
        for (int i = 0; i < grid.nN; i++) {

            sumaCP[i] = 0.0;

            for (int j = 0; j < grid.nN; j++) {
                sumaCP[i] += (hgr.C_mat[i][j] / global_data.delta_t * t0[j]);

            }

            sumaCP[i] += hgr.P_vec[i];
        }


        t0 = GaussElimination(sumaHC, sumaCP, grid.nN);

        tempMax = maxFunction(t0, grid.nN);
        tempMin = minFunction(t0, grid.nN);
        int time = step + global_data.delta_t;

       cout << "Czas: " << time << endl;
       // cout << tempMin << endl;
       //cout << "t1: " << t0 << endl;
       //for (int i = 0; i < grid.nN; i++) {
       //    cout << t0[i] << endl;
       //}
       cout<<" Temperatura minimalna: " << tempMin<< endl;
       cout << " Temperatura maksymalna: " << tempMax<< endl;

        //cout << tempMin << endl;
        //cout << tempMax << endl;
      /* zapis4 << "Czas: " << time << endl;
       zapis4 << " Temperatura minimalna: " << tempMin << endl;
       zapis4 << " Temperatura maksymalna: " << tempMax << endl;*/
       
       
    }

   // zapis4.close();
    
}
