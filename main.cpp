#if defined(UNICODE) && !defined(_UNICODE)
#define _UNICODE
#elif defined(_UNICODE) && !defined(UNICODE)
#define UNICODE
#endif

#include <tchar.h>
#include <windows.h>
#include <bits/stdc++.h>

#define ll long long
#define MAX_LOADSTRING 100
HINSTANCE hInst;								// current instance
TCHAR szTitle[MAX_LOADSTRING];					// The title bar text
TCHAR szWindowClass[MAX_LOADSTRING];			// the main window class name
COLORREF bkgndcolor;
using namespace std;

string tos(ll n)
{
    stringstream ss;
    string ans;
    ss << n;
    ss >> ans;
    return ans;
}
ll toll(string n)
{
    return atoll(n.c_str());
}

/*  Declare Windows procedure  */
LRESULT CALLBACK WindowProcedure (HWND, UINT, WPARAM, LPARAM);

/*  Make the class name into a global variable  */
TCHAR szClassName[ ] = _T("CodeBlocksWindowsApp");

int WINAPI WinMain (HINSTANCE hThisInstance,
                    HINSTANCE hPrevInstance,
                    LPSTR lpszArgument,
                    int nCmdShow)
{
    HWND hwnd;               /* This is the handle for our window */
    MSG messages;            /* Here messages to the application are saved */
    WNDCLASSEX wincl;        /* Data structure for the windowclass */

    /* The Window structure */
    wincl.hInstance = hThisInstance;
    wincl.lpszClassName = szClassName;
    wincl.lpfnWndProc = WindowProcedure;      /* This function is called by windows */
    wincl.style = CS_DBLCLKS;                 /* Catch double-clicks */
    wincl.cbSize = sizeof (WNDCLASSEX);

    /* Use default icon and mouse-pointer */
    wincl.hIcon = LoadIcon (NULL, IDI_APPLICATION);
    wincl.hIconSm = LoadIcon (NULL, IDI_APPLICATION);
    wincl.hCursor = LoadCursor (NULL, IDC_ARROW);
    wincl.lpszMenuName = NULL;                 /* No menu */
    wincl.cbClsExtra = 0;                      /* No extra bytes after the window class */
    wincl.cbWndExtra = 0;                      /* structure or the window instance */
    /* Use Windows's default colour as the background of the window */
    wincl.hbrBackground = (HBRUSH) COLOR_BACKGROUND;

    /* Register the window class, and if it fails quit the program */
    if (!RegisterClassEx (&wincl))
        return 0;

    /* The class is registered, let's create the program*/
    hwnd = CreateWindowEx (
               0,                   /* Extended possibilites for variation */
               szClassName,         /* Classname */
               _T("Code::Blocks Template Windows App"),       /* Title Text */
               WS_OVERLAPPEDWINDOW, /* default window */
               CW_USEDEFAULT,       /* Windows decides the position */
               CW_USEDEFAULT,       /* where the window ends up on the screen */
               1024,                 /* The programs width */
               765,                 /* and height in pixels */
               HWND_DESKTOP,        /* The window is a child-window to desktop */
               NULL,                /* No menu */
               hThisInstance,       /* Program Instance handler */
               NULL                 /* No Window Creation data */
           );

    /* Make the window visible on the screen */
    ShowWindow (hwnd, nCmdShow);

    /* Run the message loop. It will run until GetMessage() returns 0 */
    while (GetMessage (&messages, NULL, 0, 0))
    {
        /* Translate virtual-key messages into character messages */
        TranslateMessage(&messages);
        /* Send message to WindowProcedure */
        DispatchMessage(&messages);
    }

    /* The program return-value is 0 - The value that PostQuitMessage() gave */
    return messages.wParam;
}

///////////////////////////          line                   ///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

// simple line DDA

void LineSimpleDDA(HDC hdc, int xs, int ys, int xe, int ye, COLORREF color)
{

    double slope = (double) (ye-ys)/(xe-xs);

    if(abs(slope) <= 1)
    {
        if(xs > xe)
        {
            swap(xs,xe);
            swap(ys,ye);
        }
        int x = xs;
        int xinc = (xe-xs) > 0 ? 1:-1;
        double y = ys;
        double yinc = slope * xinc;

        while(x <= xe)
        {
            SetPixel(hdc,x,round(y),color);
            x += xinc;
            y += yinc;
        }
    }
    else
    {
        if(ys > ye)
        {
            swap(xs,xe);
            swap(ys,ye);
        }

        double iSlope = (double) 1.0/slope;
        int y = ys;
        int yinc = (ye - ys) > 0 ? 1:-1;
        double x = xs;
        double xinc = iSlope * yinc;

        while(y <= ye)
        {
            SetPixel(hdc,round(x),y,color);
            x += xinc;
            y += yinc;
        }
    }
}


void LineBresenham(HDC hdc, int xs, int ys, int xe, int ye, COLORREF color)
{
    int dx = xe - xs;
    int dy = ye - ys;

    double slope = (double) dy/dx;

    if(abs(slope) <= 1)
    {
        if(xs > xe)
        {
            std::swap(xs,xe);
            std::swap(ys,ye);
            dx *= -1;
            dy *= -1;
        }

        int x = xs, y = ys;
        int yinc = dy > 0 ? 1:-1;
        int d = dx - 2*abs(dy);
        int d1 = 2 * (dx - abs(dy));
        int d2 = -2 * abs(dy);
        while (x <= xe)
        {
            SetPixel(hdc,x,y,color);
            if(d < 0)
            {
                d += d1;
                y += yinc;
            }
            else
            {
                d += d2;
            }
            x++;
        }

    }
    else
    {
        if(ys > ye)
        {
            std::swap(xs,xe);
            std::swap(ys,ye);
            dx *= -1;
            dy *= -1;
        }

        int x = xs, y = ys;
        int xinc = dx > 0 ? 1:-1;
        int d = 2*abs(dx) - dy;
        int d1 = 2*(abs(dx) - dy);
        int d2 = 2*abs(dx);

        while(y <= ye)
        {
            SetPixel(hdc,x,y,color);
            if(d > 0)
            {
                d += d1;
                x += xinc;
            }
            else
            {
                d += d2;
            }
            y++;
        }
    }
}


void ParametricDrawLine(HDC hdc, int xs, int ys, int xe, int ye, COLORREF color)
{
    int deltax = xe - xs;
    int deltay = ye - ys;
    int n = max(abs(deltax), abs(deltay));

    double dt = 1.0 / n;
    double dx = dt * (double)deltax;
    double dy = dt * (double)deltay;

    double x = xs;
    double y = ys;

    for (int i = 0; i < n; i++)
    {
        SetPixel(hdc, round(x), round(y), color);
        x += dx;
        y += dy;
    }
}

////////////                          circle                          /////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

long dist(long x1, long y1, long x2, long y2)
{
    return sqrt(((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)));
}


void Draw8Point(HDC hdc, int xc, int yc, int a, int b, COLORREF color)
{
    SetPixel(hdc, xc + a, yc + b, color);
    SetPixel(hdc, xc - a, yc + b, color);
    SetPixel(hdc, xc + a, yc - b, color);
    SetPixel(hdc, xc - a, yc - b, color);
    SetPixel(hdc, xc + b, yc + a, color);
    SetPixel(hdc, xc - b, yc + a, color);
    SetPixel(hdc, xc + b, yc - a, color);
    SetPixel(hdc, xc - b, yc - a, color);
}


void CircleDirectCartisian(HDC hdc, int xc, int yc, int R, COLORREF color)
{
    int x = 0, y = R;
    int R2 = R*R;

    while(x <= y)
    {
        Draw8Point(hdc,xc,yc,x,y,color);
        x++;
        y = round(sqrt((double) R2 - x*x));
    }
}




void CircleDirectPolar(HDC hdc, int xc, int yc, int R, COLORREF color)
{
    int x = R, y = 0;
    double theta = 0, dtheta = 1.0/R;

    while(x >= y)
    {
        Draw8Point(hdc,xc,yc,x,y,color);
        theta += dtheta;
        x = round(R * cos(theta));
        y = round(R * sin(theta));
    }
}



void CircleIterativePolar(HDC hdc, int xc, int yc, int R, COLORREF color)
{
    double x = R, y = 0;
    double dtheta = 1.0/R;
    double cos_dtheta = cos(dtheta);
    double sin_dtheta = sin(dtheta);

    while(x >= y)
    {
        Draw8Point(hdc,xc,yc,round(x),round(y),color);
        x = x*cos_dtheta - y*sin_dtheta;
        y = x*sin_dtheta + y*cos_dtheta;
    }
}


void CircleBresenham(HDC hdc, int xc, int yc, int R, COLORREF color)
{
    int x = 0, y = R;
    int d = 1 - R;

    while(x <= y)
    {
        Draw8Point(hdc,xc,yc,round(x),round(y),color);
        if(d < 0)
        {
            d += 2*x + 3; // (x+1,y) ... next midpoint is inside the circle
        }
        else
        {
            d += 2*(x-y)+5;  //(x+1,y-1) ... next midpoint is outside the circle
            y--;
        }
        x++;
    }
}


void CircleFasterBresenham(HDC hdc, int xc, int yc, int R, COLORREF color)
{
    int x = 0, y = R;
    int d = 1-R;
    int c1 = 3;
    int c2 = 5 - 2*R;

    while(x <= y)
    {
        Draw8Point(hdc,xc,yc,round(x),round(y),color);
        if(d < 0)
        {
            d += c1; // (x+1,y) ... next midpoint is inside the circle
            c2 += 2;
        }
        else
        {
            d += c2;  //(x+1,y-1) ... next midpoint is outside the circle
            c2 += 4;
            y--;
        }
        c1 += 2;
        x++;
    }
}



////////////                          curve                           /////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

struct Vector2
{

    double x,y;
    Vector2(double a = 0, double b = 0)
    {
        x = a;
        y = b;
    }

};

class Vector4
{
    double v[4]= {0.0,0.0,0.0,0.0};
public:

    Vector4(double a = 0, double b = 0, double c = 0, double d = 0)
    {
        v[0] = a;
        v[1] = b;
        v[2] = c;
        v[3] = d;
    }

    Vector4(double a[])
    {
        memcpy(v,a, 4 * sizeof(double));
    }

    double & operator[](int i)
    {
        return v[i];
    }
};


class Matrix4
{
    Vector4 M[4];
public:

    Matrix4(double A[])
    {
        memcpy(M,A, 16 * sizeof(double));
    }

    Vector4 & operator[](int i)
    {
        return M[i];
    }
};


Vector4 & operator *(Matrix4 M, Vector4 &b)
{
    Vector4 res;
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4 ; j++)
        {
            res[i] += M[i][j] * b[j];
        }
    }

    return res;
}


double dotproduct(Vector4 & a, Vector4 & b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
}


Vector4 getHermitecoeff(double x0, double s0, double x1, double s1)
{
    static double H[16] = {2, 1, -2, 1, -3, -2, 3, -1, 0, 1, 0, 0, 1, 0, 0, 0};
    static Matrix4 basisMatrix(H);
    Vector4 V(x0,s0,x1,s1);

    return basisMatrix * V;
}


void DrawHermiteCurve(HDC hdc, Vector2 &p0, Vector2 &T0, Vector2 &p1, Vector2 &T1, int numofpoints, COLORREF color)
{
    Vector4 xCoeff = getHermitecoeff(p0.x, T0.x, p1.x, T1.x);
    Vector4 yCoeff = getHermitecoeff(p0.y, T0.y, p1.y, T1.y);

    if(numofpoints < 2) return;
    double dt = 1.0/(numofpoints-1);
    for(double t = 0 ; t <= 1 ; t += dt)
    {
        Vector4 vt;
        vt[3] = 1;
        for(int i = 2 ; i >= 0 ; i--)
        {
            vt[i] = vt[i+1] * t;
        }


        int x = round(dotproduct(xCoeff, vt));
        int y = round(dotproduct(yCoeff, vt));

        if(t == 0) MoveToEx(hdc,x,y,NULL);
        else LineTo(hdc,x,y);
    }

}


void DrawBezierCurve(HDC hdc,Vector2& P0,Vector2& P1,Vector2& P2,Vector2& P3,int numpoints , COLORREF c)
{
    Vector2 T0(3*(P1.x-P0.x),3*(P1.y-P0.y));
    Vector2 T1(3*(P3.x-P2.x),3*(P3.y-P2.y));
    DrawHermiteCurve(hdc,P0,T0,P3,T1,numpoints, c);
}


void DrawCardinalSpline(HDC hdc,Vector2 P[],int n,int numpix, COLORREF co)
{
    double c1=0.5;
    Vector2 T0(c1*(P[2].x-P[0].x),c1*(P[2].y-P[0].y));
    for(int i=2; i<n-1; i++)
    {
        Vector2 T1(c1*(P[i+1].x-P[i-1].x),c1*(P[i+1].y-P[i-1].y));
        DrawHermiteCurve(hdc,P[i-1],T0,P[i],T1,numpix , co);
        T0=T1;
    }
}

////////////                          filling                         /////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
#define MAXENTRIES 600


struct Entry
{
    int xmin, xmax;
};

void initEntries(Entry table[])
{
    for(int i = 0 ; i < MAXENTRIES ; i++)
    {
        table[i].xmin = INT_MAX;
        table[i].xmax = INT_MIN;
    }
}

void scanEdge(POINT v1, POINT v2, Entry table[])
{
    if(v1.y == v2.y) return;
    if(v1.y > v2.y)
    {
        std::swap(v1,v2);
    }

    double invSlope = (double)(v2.x - v1.x)/(v2.y - v1.y);
    double x = v1.x;
    int y = v1.y;
    while(y < v2.y)
    {
        if(x < table[y].xmin) table[y].xmin = (int)ceil(x);
        if(x > table[y].xmax) table[y].xmax = (int)floor(x);
        y++;
        x += invSlope;
    }
}

void drawScanLines(HDC hdc, Entry table[], COLORREF color)
{
    for(int y = 0 ; y < MAXENTRIES ; y++)
    {
        if(table[y].xmin < table[y].xmax)
        {
            for(int x = table[y].xmin ; x <= table[y].xmax ; x++)
            {
                SetPixel(hdc,x,y,color);
            }
        }
    }
}

void convexFill(HDC hdc, POINT p[], int n, COLORREF color)
{
    Entry *table = new Entry[MAXENTRIES];
    initEntries(table);
    POINT v1 = p[n-1];
    for(int i = 0 ; i < n ; i++)
    {
        POINT v2 = p[i];
        scanEdge(v1,v2,table);
        v1 = p[i];
    }
    drawScanLines(hdc,table,color);
    delete table;

}

/////////////// clipping//////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

union OutCode
{
    unsigned All:4;
    struct
    {
        unsigned left:1,top:1,right:1,bottom:1;
    };
};

OutCode GetOutCode(double x, double y, int xleft,int ytop,int xright,int ybottom)
{
    OutCode out;
    out.All=0;
    if(x < xleft) out.left = 1;
    else if(x > xright) out.right = 1;
    if(y < ytop) out.top = 1;
    else if(y > ybottom) out.bottom = 1;
    return out;
}

void VIntersect(double xs,double ys,double xe,double ye,int x,double *xi,double *yi)
{
    *xi = x; // x stands for xEdge .. left or right edges
    *yi = ys+(x-xs)*(ye-ys)/(xe-xs);
}
void HIntersect(double xs,double ys,double xe,double ye,int y,double *xi,double *yi)
{
    *yi=y;
    *xi=xs+(y-ys)*(xe-xs)/(ye-ys);
}

void lineClipping(HDC hdc,int xs,int ys,int xe,int ye,int xleft,int ytop,int xright,int ybottom)
{
    double x1 = xs, y1 = ys, x2 = xe, y2 = ye;
    OutCode out1 = GetOutCode(x1,y1, xleft,ytop,xright,ybottom);
    OutCode out2=GetOutCode(x2,y2,xleft,ytop,xright,ybottom);
    //0000
    //0000


    while((out1.All || out2.All) && !(out1.All & out2.All) ) // check condition
    {
        double xi, yi;
        if(out1.All != 0)
        {
            if(out1.left) VIntersect(x1,y1,x1,y2,xleft,&xi, &yi);
            else if(out1.top)HIntersect(x1,y1,x2,y2,ytop,&xi,&yi);
            else if(out1.right)VIntersect(x1,y1,x2,y2,xright,&xi,&yi);
            else HIntersect(x1,y1,x2,y2,ybottom,&xi,&yi);
            x1 = xi;
            y1 = yi;
            out1=GetOutCode(x1,y1,xleft,ytop,xright,ybottom);
        }
        else
        {
            if(out2.left)VIntersect(x1,y1,x2,y2,xleft,&xi,&yi);
            else if(out2.top)HIntersect(x1,y1,x2,y2,ytop,&xi,&yi);
            else if(out2.right)VIntersect(x1,y1,x2,y2,xright,&xi,&yi);
            else HIntersect(x1,y1,x2,y2,ybottom,&xi,&yi);
            x2 = xi;
            y2 = yi;
            out2 = GetOutCode(x2,y2,xleft,ytop,xright,ybottom);
        }
    }
    if(!out1.All && !out2.All)
    {
        MoveToEx(hdc,round(x1),round(y1),NULL);
        LineTo(hdc,round(x2),round(y2));
    }
}

/////////////////////////////////////////////////////
/////////////////Bounce ////////////////////////////

void circlePointClipping(HDC hdc, int x , int y, int xc, int yc, int R, COLORREF color)
{
    double distance = sqrt( ( (x-xc)*(x-xc) ) + ( (y-yc)*(y-yc) ) );

    if(distance < R)
    {
        CircleDirectPolar(hdc,xc,yc,R,color);
        SetPixel(hdc, x,y,color);
    }
}


void CircleLineClipping(HDC hdc, int xs, int ys, int xe, int ye, int xc, int yc, int R, COLORREF color)
{

    double slope = (double) (ye-ys)/(xe-xs);
    CircleDirectPolar(hdc,xc,yc,R,color);
    if(abs(slope) <= 1)
    {
        if(xs > xe)
        {
            swap(xs,xe);
            swap(ys,ye);
        }
        int x = xs;
        int xinc = (xe-xs) > 0 ? 1:-1;
        double y = ys;
        double yinc = slope * xinc;

        while(x <= xe)
        {
            int ytmp = round(y);
            double distance = sqrt( ( (x-xc)*(x-xc) ) + ( (ytmp-yc)*(ytmp-yc) ) );
            if(distance < R)
            {
                SetPixel(hdc, x, ytmp, color);
            }
            x += xinc;
            y += yinc;
        }
    }
    else
    {
        if(ys > ye)
        {
            swap(xs,xe);
            swap(ys,ye);
        }

        double iSlope = (double) 1.0/slope;
        int y = ys;
        int yinc = (ye - ys) > 0 ? 1:-1;
        double x = xs;
        double xinc = iSlope * yinc;

        while(y <= ye)
        {
            int xtmp = round(x);
            double distance = sqrt( ( (xtmp-xc)*(xtmp-xc) ) + ( (y-yc)*(y-yc) ) );
            if(distance < R)
            {
                SetPixel(hdc, xtmp, y, color);
            }
            x += xinc;
            y += yinc;
        }
    }
}




void PointClipping(HDC hdc,int x,int y,int xleft,int ytop,int xright,int ybottom,COLORREF color)
{
if(x>=xleft && x<= xright && y>=ytop && y<=ybottom)
    SetPixel(hdc,x,y,color);
}






//////////////////////////////////////////////////////////////////////////

#define ID_Load 100
#define ID_RED 0
#define ID_GREEN 1
#define ID_BLUE 2
#define ID_Exit 3
#define ID_Save 4
#define ID_DDA 5
#define ID_MidpointLine 6
#define ID_ParametricLine 7
#define ID_Cartesian 8
#define ID_Polar 9
#define ID_ITerativePolar 10
#define ID_ParametricCircle 11
#define ID_MidPointCircle 12
#define ID_FastMidPoint 13

#define ID_FirstDegree 14
#define ID_SecondDegree 15
#define ID_Hermite 16
#define ID_Bezier 17
#define ID_Spline 18

#define ID_ConvexFilling 19

#define ID_ClippingLine 20
#define ID_ClippingPoint 21
#define ID_ClippingCircle 22
#define ID_ClippingRectanlePoint 23


/*  This function is called by the Windows function DispatchMessage()  */


void ChangeBkGnd(HWND hWnd,int ch)
{
    if(ch == 0)bkgndcolor=RGB(255,0,0);
    else if(ch == 1)bkgndcolor=RGB(0,255,0);
    else if(ch == 2)bkgndcolor=RGB(0,0,255);
    InvalidateRect(hWnd,NULL, TRUE);
}

string LineEncode(int ID , int xs,int ys,int xe, int ye)
{
    return tos(ID) + '/' + tos(xs) + '/' + tos(ys) + '/' + tos(xe) + '/' + tos(ye) + '/';
}

string CircleEncode (int ID , int x,int y,int radius)
{
    return tos(ID)+'/'+tos(x)+'/'+tos(y)+'/'+tos(radius)+'/';
}


string CurveFourPoint (int ID , POINT p1,POINT p2,POINT p3,POINT p4)
{
    return tos(ID)+'/'+tos(p1.x)+'/'+tos(p1.y)+'/'+tos(p2.x)+'/'+tos(p2.y)+
            '/'+tos(p3.x)+'/'+tos(p3.y)+'/'+tos(p4.x)+'/'+tos(p4.y)+'/';
}
string CurveFour (int ID , Vector2 p1,Vector2 p2,Vector2 p3,Vector2 p4)
{
    return tos(ID)+'/'+tos(p1.x)+'/'+tos(p1.y)+'/'+tos(p2.x)+'/'+tos(p2.y)+
            '/'+tos(p3.x)+'/'+tos(p3.y)+'/'+tos(p4.x)+'/'+tos(p4.y)+'/';
}

string CurveSix (int ID , Vector2 p1,Vector2 p2,Vector2 p3,Vector2 p4,Vector2 p5,Vector2 p6)
{
    return tos(ID)+'/'+tos(p1.x)+'/'+tos(p1.y)+'/'+tos(p2.x)+'/'+tos(p2.y)
            +'/'+tos(p3.x)+'/'+tos(p3.y)+'/'+tos(p4.x)+'/'+tos(p4.y)
            +'/'+tos(p5.x)+'/'+tos(p5.y)+'/'+tos(p6.x)+'/'+tos(p6.y)+'/';
}

string CurveSixPoint (int ID , POINT p1,POINT p2,POINT p3,POINT p4,POINT p5,POINT p6)
{
    return tos(ID)+'/'+tos(p1.x)+'/'+tos(p1.y)+'/'+tos(p2.x)+'/'+tos(p2.y)
            +'/'+tos(p3.x)+'/'+tos(p3.y)+'/'+tos(p4.x)+'/'+tos(p4.y)
            +'/'+tos(p5.x)+'/'+tos(p5.y)+'/'+tos(p6.x)+'/'+tos(p6.y)+'/';
}

string CurveThreePoint (int ID , POINT p1,POINT p2,POINT p3)
{
    return tos(ID)+'/'+tos(p1.x)+'/'+tos(p1.y)+'/'+tos(p2.x)+'/'+tos(p2.y)
            +'/'+tos(p3.x)+'/'+tos(p3.y)+'/';
}


vector<int> enCode(string x)
{
    string tmp = "";
    vector<int> integars;
    for(int i=0; i<x.size(); i++)
    {
        if(x[i] == '/')
        {
            integars.push_back(toll(tmp));
            tmp = "";
        }
        else tmp+=x[i];
    }
    return integars;
}

void DoLoad(vector<string> file , HDC hdc)
{

    for(int i = 0 ; i <file.size(); i++)
    {
        vector<int> code = enCode(file[i]);
        if(code[0] == ID_Cartesian)
        {
            CircleDirectCartisian(hdc, code[1], code[2], code[3],RGB(255,255,255));
        }
        else if(code[0] == ID_Polar)
        {
            CircleDirectPolar(hdc, code[1], code[2], code[3], RGB(255,255,255));
        }
        else if(code[0] == ID_ITerativePolar)
        {
            CircleIterativePolar(hdc, code[1], code[2], code[3],RGB(255,255,255));
        }
        else if(code[0] == ID_MidPointCircle)
        {
            CircleBresenham(hdc, code[1], code[2], code[3], RGB(255,255,255));
        }
        else if(code[0] == ID_FastMidPoint)
        {
            CircleFasterBresenham(hdc,  code[1], code[2], code[3], RGB(255,255,255));
        }
        //line
        else if(code[0] == ID_DDA)
        {
            LineSimpleDDA(hdc,code[1], code[2], code[3],code[4],RGB(255,255,255));
        }
        else if(code[0] == ID_MidpointLine)
        {
            LineBresenham(hdc,code[1], code[2], code[3],code[4],RGB(255,255,255));
        }
        else if(code[0] == ID_ParametricLine)
        {
            ParametricDrawLine(hdc,code[1], code[2], code[3],code[4],RGB(255,255,255));
        }
        else if(code[0] == ID_SecondDegree || code[0] ==ID_Bezier){
            Vector2 p[4];
            int cnt=0;
            for(int i=1;code.size();i+=2){
                p[cnt].x = code[i];
                p[cnt].y = code[i+1];
                cnt++;
            }
            DrawBezierCurve(hdc, p[0],p[1],p[2],p[3], 50, RGB(255,255,255));
        }
        else if(code[0] == ID_Hermite){
            Vector2 p[4];
            int cnt=0;
            for(int i=1;code.size();i+=2){
                p[cnt].x = code[i];
                p[cnt].y = code[i+1];
                cnt++;
            }
            DrawHermiteCurve(hdc, p[0],p[1],p[2],p[3], 50, RGB(255,255,255));
        }else if(code[0] == ID_Spline){
                Vector2 p[6];
                int indx = 0;
                for(int i=1;i<code.size();i+=2){
                    p[indx].x = code[i];
                    p[indx].y = code[i+1];
                    indx++;
                }
                DrawCardinalSpline(hdc, p,6, 50, RGB(255,255,255));
        }else if(code[0] == ID_ClippingPoint){
            POINT p[3];
            int ind = 0;
            for(int i=1;i<code.size();i+=2){
                p[ind].x = code[i];
                p[ind].y = code[i+1];
                ind++;
            }
            circlePointClipping(hdc,p[2].x,p[2].y,p[0].x,p[0].y,dist(p[0].x,p[0].y,p[1].x,p[1].y),RGB(255,255,255));
        }else if(code[0] == ID_ClippingLine){
            POINT points[4];
            int ind = 0;
            for(int i=1;i<code.size();i+=2){
                points[ind].x = code[i];
                points[ind].y = code[i+1];
                ind++;
            }
            Rectangle(hdc, points[0].x, points[1].y, points[1].x,points[0].y);
            lineClipping(hdc,points[2].x,points[2].y,points[3].x,points[3].y,points[0].x, points[1].y, points[1].x,points[0].y);
        }
    }
}



LRESULT CALLBACK WindowProcedure (HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam)
{

    static int flag = 0;
    static vector<POINT> points;
    int wmId, wmEvent;
    static PAINTSTRUCT ps;
    static HDC hdc;
    static HBRUSH hBrush = CreateSolidBrush(RGB(230,230,230));
    static int choice ;
    static vector<string> file;

    switch (message)                  /* handle the messages */
    {
    case WM_LBUTTONDOWN:
            POINT tmp;
            tmp.x = LOWORD(lParam);
            tmp.y = HIWORD(lParam);
            flag++;
            points.push_back(tmp);
            break;
        break;
    case WM_CREATE:
    {
        HMENU hMenubar = CreateMenu();
        HMENU hLine = CreateMenu();
        HMENU hCircle = CreateMenu();
        HMENU hCurve = CreateMenu();
        HMENU hColor = CreateMenu();
        HMENU hFile  = CreateMenu();
        HMENU hFilling = CreateMenu();
        HMENU hClipping = CreateMenu();

        AppendMenu(hMenubar, MF_POPUP, (UINT_PTR)hFile, "File");
        AppendMenu(hMenubar, MF_POPUP, (UINT_PTR)hLine, "Line");
        AppendMenu(hMenubar, MF_POPUP, (UINT_PTR)hCircle, "Circle");
        AppendMenu(hMenubar, MF_POPUP, (UINT_PTR)hCurve, "Curve");
        AppendMenu(hMenubar, MF_POPUP, (UINT_PTR)hClipping, "Clipping");
        AppendMenu(hMenubar, MF_POPUP, (UINT_PTR)hFilling, "Convex filling");
        AppendMenu(hMenubar, MF_POPUP, (UINT_PTR)hColor, "Color BackGround");

        AppendMenu(hFilling, MF_STRING, ID_ConvexFilling, "Convex filling");

        AppendMenu(hClipping, MF_STRING, ID_ClippingLine, "Clipping rectangle line");
        AppendMenu(hClipping, MF_STRING, ID_ClippingRectanlePoint, "Clipping rectangle point");
        AppendMenu(hClipping, MF_STRING, ID_ClippingCircle, "Clipping Circle line");
        AppendMenu(hClipping, MF_STRING, ID_ClippingPoint, "Clipping Circle Point");

        AppendMenu(hLine, MF_STRING, ID_DDA, "DDA");
        AppendMenu(hLine, MF_STRING, ID_MidpointLine, "Midpoint");
        AppendMenu(hLine, MF_STRING, ID_ParametricLine, "Parametric");


        AppendMenu(hCircle, MF_STRING, ID_Cartesian, "Cartesian");
        AppendMenu(hCircle, MF_STRING, ID_Polar, "Polar");
        AppendMenu(hCircle, MF_STRING, ID_ITerativePolar, "Iterative Polar");
        AppendMenu(hCircle, MF_STRING, ID_MidPointCircle, "Midpoint");
        AppendMenu(hCircle, MF_STRING, ID_FastMidPoint, "Fast Midpoint");
        //AppendMenu(hCircle, MF_STRING, ID_Polar, "Parametric Circle");

        AppendMenu(hColor, MF_STRING, ID_RED, "Red");
        AppendMenu(hColor, MF_STRING, ID_GREEN, "Green");
        AppendMenu(hColor, MF_STRING, ID_BLUE, "Blue");

        AppendMenu(hCurve, MF_STRING, ID_FirstDegree, "First Degree");
        AppendMenu(hCurve, MF_STRING, ID_SecondDegree, "Second Degree");
        AppendMenu(hCurve, MF_STRING, ID_Hermite, "Third Hermite");
        AppendMenu(hCurve, MF_STRING, ID_Bezier, "Third Bezier");
        AppendMenu(hCurve, MF_STRING, ID_Spline, "Third Spline");

        AppendMenu(hFile, MF_STRING, ID_Save, "Save");
        AppendMenu(hFile, MF_STRING, ID_Load, "Load");
        AppendMenu(hFile, MF_STRING, ID_Exit, "Exit");

        SetMenu(hwnd, hMenubar);
        break;
    }

    case WM_COMMAND:
    {
        wmId    = LOWORD(wParam);
        wmEvent = HIWORD(wParam);
        if(LOWORD(wParam) == ID_Exit)
        {
            exit(0);
            break;
        }
        else if (LOWORD(wParam) == ID_RED)
        {
            ChangeBkGnd(hwnd,0);
        }
        else if(LOWORD(wParam) == ID_GREEN)
        {
            ChangeBkGnd(hwnd,1);
        }
        else if(LOWORD(wParam) == ID_BLUE)
        {
            ChangeBkGnd(hwnd,2);
        }
        else
        {
            choice = LOWORD(wParam);
        }
    }
    case WM_ERASEBKGND:
        HPEN pen;
        HBRUSH	brush;
        RECT rect;

        pen=CreatePen(PS_SOLID, 1, bkgndcolor);
        brush=CreateSolidBrush(bkgndcolor);
        SelectObject((HDC)wParam, pen);
        SelectObject((HDC)wParam, brush);

        GetClientRect(hwnd, &rect);

        Rectangle((HDC)wParam, rect.left, rect.top, rect.right, rect.bottom);

        break;
    case WM_PAINT:
        hdc = GetDC(hwnd);

        if(choice == ID_Save)
        {
            cout<<"Save";
            ofstream out("saveData.txt");
            out<<file.size()<<"\n";
            for(int i=0; i<file.size(); i++)
            {
                out << file[i]<<"\n";
            }
            out.close();
            choice = -9999999;
        }
        if(choice == ID_Load)
        {
            cout<<"LOADING";
            ifstream fin("saveData.txt");
            if(fin.is_open())
            {
                int sz;
                fin>>sz;
                string tmp;
                for(int i = 0; i < sz; ++i)
                {
                    fin >> tmp;
                    file.push_back(tmp);
                    cout<<file[i]<<"\n";
                }
            }
            choice = -9999999;
            DoLoad(file,hdc);
        }

        if(flag == 2)
        {
            cout<<choice;
            if(choice == ID_Cartesian)
            {
                CircleDirectCartisian(hdc, points[0].x,points[0].y , dist(points[0].x,points[0].y,points[1].x,points[1].y),RGB(255,255,255));
                file.push_back(CircleEncode(ID_Cartesian , points[0].x,points[0].y , dist(points[0].x,points[0].y,points[1].x,points[1].y)));
                flag = 0;points.clear();
            }
            else if(choice == ID_Polar)
            {
                CircleDirectPolar(hdc,  points[0].x,points[0].y , dist(points[0].x,points[0].y,points[1].x,points[1].y), RGB(255,255,255));
                file.push_back(CircleEncode(ID_Polar , points[0].x,points[0].y , dist(points[0].x,points[0].y,points[1].x,points[1].y)));
                flag = 0;points.clear();
            }
            else if(choice == ID_ITerativePolar)
            {
                CircleIterativePolar(hdc,  points[0].x,points[0].y , dist(points[0].x,points[0].y,points[1].x,points[1].y),RGB(255,255,255));
                file.push_back(CircleEncode(ID_ITerativePolar , points[0].x,points[0].y , dist(points[0].x,points[0].y,points[1].x,points[1].y)));
                flag = 0;points.clear();
            }
            else if(choice == ID_MidPointCircle)
            {
                CircleBresenham(hdc, points[0].x,points[0].y , dist(points[0].x,points[0].y,points[1].x,points[1].y), RGB(255,255,255));
                file.push_back(CircleEncode(ID_MidPointCircle , points[0].x,points[0].y , dist(points[0].x,points[0].y,points[1].x,points[1].y)));
                flag = 0;points.clear();
            }
            else if(choice == ID_FastMidPoint)
            {
                CircleFasterBresenham(hdc,  points[0].x,points[0].y , dist(points[0].x,points[0].y,points[1].x,points[1].y), RGB(255,255,255));
                file.push_back(CircleEncode(ID_FastMidPoint , points[0].x,points[0].y , dist(points[0].x,points[0].y,points[1].x,points[1].y)));
                flag = 0;points.clear();
            }
            else if(choice == ID_DDA)
            {
                LineSimpleDDA(hdc,points[0].x,points[0].y,points[1].x,points[1].y,RGB(255,255,255));
                file.push_back(LineEncode(ID_DDA , points[0].x,points[0].y, points[1].x,points[1].y));
                flag = 0;points.clear();
            }
            else if(choice == ID_MidpointLine)
            {
                LineBresenham(hdc,points[0].x,points[0].y,points[1].x,points[1].y,RGB(255,255,255));
                file.push_back(LineEncode(ID_MidpointLine ,points[0].x,points[0].y,points[1].x,points[1].y));
                flag = 0;points.clear();
            }
            else if (choice ==  ID_ParametricLine || choice == ID_FirstDegree)
            {
                ParametricDrawLine(hdc,points[0].x,points[0].y,points[1].x,points[1].y,RGB(255,255,255));
                file.push_back(LineEncode(ID_ParametricLine ,points[0].x,points[0].y,points[1].x,points[1].y));
                flag = 0;points.clear();
            }

        }

        if(flag == 3)
        {
            if(choice == ID_SecondDegree)
            {
                Vector2 p1;
                p1.x = points[0].x*1.0;
                p1.y = points[0].y*1.0;
                Vector2 p2;
                p2.x = points[1].x*1.0;
                p2.y = points[1].y*1.0;
                Vector2 p3;
                p3.x = points[2].x*1.0;
                p3.y = points[2].y*1.0;
                DrawBezierCurve(hdc, p1,p2,p2,p3, 50, RGB(255,255,255));
                file.push_back(CurveFour(ID_SecondDegree,p1,p2,p2,p3));
                flag=0;
                points.clear();
            }else if(choice == ID_ClippingPoint){
                circlePointClipping(hdc,points[2].x,points[2].y,points[0].x,points[0].y,dist(points[0].x,points[0].y,points[1].x,points[1].y),RGB(255,255,255));
                file.push_back(CurveThreePoint(ID_ClippingPoint ,points[0],points[1],points[2] ));
                flag=0;
                points.clear();
            }
            else if(choice == ID_ClippingRectanlePoint){
                Rectangle(hdc, points[0].x, points[1].y, points[1].x,points[0].y);
                PointClipping(hdc,points[2].x,points[2].y,points[0].x, points[1].y, points[1].x,points[0].y, RGB(255,0,255));
                file.push_back(CurveThreePoint(ID_ClippingPoint ,points[0],points[1],points[2] ));
                flag=0;
                points.clear();
            }
        }

        if(flag == 4)
        {
                Vector2 p1;
                p1.x = points[0].x*1.0;
                p1.y = points[0].y*1.0;
                Vector2 p2;
                p2.x = points[1].x*1.0;
                p2.y = points[1].y*1.0;
                Vector2 p3;
                p3.x = points[2].x*1.0;
                p3.y = points[2].y*1.0;
                Vector2 p4;
                p4.x = points[3].x*1.0;
                p4.y = points[3].y*1.0;
            if(choice == ID_Hermite)
            {

                DrawHermiteCurve(hdc, p1,p2,p3,p4, 50, RGB(255,255,255));
                file.push_back(CurveFour(ID_SecondDegree,p1,p2,p3,p4));
                flag=0;
                points.clear();
            }
            else if(choice == ID_Bezier)
            {
                DrawBezierCurve(hdc, p1,p2,p3,p4, 50, RGB(255,255,255));
                 file.push_back(CurveFour(ID_SecondDegree,p1,p2,p3,p4));
                flag=0;
                points.clear();
            }else if(choice == ID_ClippingLine){
                 //first point is the botton left and seconed point in the top right
                 //don't forget
                 Rectangle(hdc, points[0].x, points[1].y, points[1].x,points[0].y);
                 lineClipping(hdc,points[2].x,points[2].y,points[3].x,points[3].y,points[0].x, points[1].y, points[1].x,points[0].y);
                 file.push_back(CurveFourPoint(ID_ClippingLine,points[0],points[1],points[2],points[3]));
                 flag=0;
                 points.clear();
            }else if(choice == ID_ClippingCircle){
                CircleLineClipping(hdc,points[2].x,points[2].y,points[3].x,points[3].y,points[0].x,points[0].y,dist(points[0].x,points[0].y,points[1].x,points[1].y),RGB(255,255,255));
                file.push_back(CurveFourPoint(ID_ClippingLine,points[0],points[1],points[2],points[3]));
                flag=0;
                points.clear();
            }

        }

        if(flag == 6){
            if(choice == ID_Spline)
            {   Vector2 p[6];
                for(int i=0;i<6;i++){
                    p[i].x = points[i].x*1.0;
                    p[i].y = points[i].y*1.0;
                }
                DrawCardinalSpline(hdc, p,6, 50, RGB(255,255,255));
                file.push_back(CurveSix(ID_Spline ,p[0],p[1],p[2],p[3],p[4],p[5] ));
                flag=0;
                points.clear();
            }else if(choice == ID_ConvexFilling){
                POINT p[6];
                for(int i=0;i<6;i++){
                    p[i].x = points[i].x;
                    p[i].y = points[i].y;
                }
                convexFill(hdc, p,6, RGB(255,255,255));
                file.push_back(CurveSixPoint(ID_Spline ,p[0],p[1],p[2],p[3],p[4],p[5] ));
                flag=0;
                points.clear();
            }
        }

        EndPaint(hwnd, &ps);
        break;

    case WM_DESTROY:
        PostQuitMessage (0);       /* send a WM_QUIT to the message queue */
        break;
    default:                      /* for messages that we don't deal with */
        return DefWindowProc (hwnd, message, wParam, lParam);
    }

    return 0;
}
