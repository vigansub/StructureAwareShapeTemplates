#include "matlabwrapper.h"
#include "engine.h"
#include <QDebug>
#define BUFSIZE 256

using namespace std;

void matlab::example_1()
{
    Engine *me = engOpen(NULL);
    if (!me)
    {
        qDebug() << "Can't start Matlab engine!";
    }
    engSetVisible(me, true);

    const int N= 100;
    vector<double> x(N);
    iota(x.begin(), x.end(), 0);
    transform(x.begin(), x.end(), x.begin(), [](double d){return d/10.;});

    mxArray *mx_x = mxCreateDoubleMatrix(1, N, mxREAL);
    memcpy(mxGetPr(mx_x), &x.front(), N*sizeof(double));

    engPutVariable(me, "x", mx_x);
    engEvalString(me, "y = sin(x);");
    engEvalString(me, "plot(x, y);");

    fgetc(stdin);

    mxDestroyArray(mx_x);
    engClose(me);
}

void matlab::example_2()
{
    Engine *ep;
    mxArray *T = NULL, *result = NULL;
    char buffer[BUFSIZE+1];
    double time[10] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };

    /*
     * Call engOpen with a NULL string. This starts a MATLAB process
     * on the current host using the command "matlab".
     */
    if (!(ep = engOpen(""))) {
        fprintf(stderr, "\nCan't start MATLAB engine\n");
        return;
    }

    /*
     * PART I
     *
     * For the first half of this demonstration, we will send data
     * to MATLAB, analyze the data, and plot the result.
     */

    /*
     * Create a variable for our data
     */
    T = mxCreateDoubleMatrix(1, 10, mxREAL);
    memcpy((void *)mxGetPr(T), (void *)time, sizeof(time));
    /*
     * Place the variable T into the MATLAB workspace
     */
    engPutVariable(ep, "T", T);

    /*
     * Evaluate a function of time, distance = (1/2)g.*t.^2
     * (g is the acceleration due to gravity)
     */
    engEvalString(ep, "D = .5.*(-9.8).*T.^2;");

    /*
     * Plot the result
     */
    engEvalString(ep, "plot(T,D);");
    engEvalString(ep, "title('Position vs. Time for a falling object');");
    engEvalString(ep, "xlabel('Time (seconds)');");
    engEvalString(ep, "ylabel('Position (meters)');");

    /*
     * use fgetc() to make sure that we pause long enough to be
     * able to see the plot
     */
    printf("Hit return to continue\n\n");
    fgetc(stdin);
    /*
     * We're done for Part I! Free memory, close MATLAB figure.
     */
    printf("Done for Part I.\n");
    mxDestroyArray(T);
    engEvalString(ep, "close;");


    /*
     * PART II
     *
     * For the second half of this demonstration, we will request
     * a MATLAB string, which should define a variable X.  MATLAB
     * will evaluate the string and create the variable.  We
     * will then recover the variable, and determine its type.
     */

    /*
     * Use engOutputBuffer to capture MATLAB output, so we can
     * echo it back.  Ensure first that the buffer is always NULL
     * terminated.
     */

    buffer[BUFSIZE] = '\0';
    engOutputBuffer(ep, buffer, BUFSIZE);
    while (result == NULL) {
        char str[BUFSIZE+1];
        /*
         * Get a string input from the user
         */
        printf("Enter a MATLAB command to evaluate.  This command should\n");
        printf("create a variable X.  This program will then determine\n");
        printf("what kind of variable you created.\n");
        printf("For example: X = 1:5\n");
        printf(">> ");

        fgets(str, BUFSIZE, stdin);

        /*
         * Evaluate input with engEvalString
         */
        engEvalString(ep, str);

        /*
         * Echo the output from the command.
         */
        printf("%s", buffer);

        /*
         * Get result of computation
         */
        printf("\nRetrieving X...\n");
        if ((result = engGetVariable(ep,"X")) == NULL)
          printf("Oops! You didn't create a variable X.\n\n");
        else {
        printf("X is class %s\t\n", mxGetClassName(result));
        }
    }

    /*
     * We're done! Free memory, close MATLAB engine and exit.
     */
    printf("Done!\n");
    mxDestroyArray(result);
    engClose(ep);
}

// to test different built in functions of Matlab
void matlab::example_3()
{
    Engine *ep;
    if (!(ep = engOpen(""))) {
        fprintf(stderr, "\nCan't start MATLAB engine\n");
        return;
    }
    engEvalString(ep, "x=1:0.1:10;");
    engEvalString(ep, "y=sin(x);");
    engEvalString(ep, "plot(x,y);");
    engEvalString(ep, "dejong5fcn;");
    engClose(ep);
}

// test simple optimization procedure
void matlab::example_4()
{
    Engine *ep;
    if (!(ep = engOpen(""))) {
        fprintf(stderr, "\nCan't start MATLAB engine\n");
        return;
    }
    engEvalString(ep, "fun = @(x)100*(x(2)-x(1)^2)^2 + (1-x(1))^2;"
                      "x0 = [-1,2];"
                      "A = [1,2];"
                      "b = 1;"
                      "x = fmincon(fun,x0,A,b)");
    engClose(ep);
}
