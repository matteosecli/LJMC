/* pBar - An extremely easy way to get a higly customizable
 * progress bar in a C++ command line application
 * (Mac, Windows, & Linux).
 *
 * Copyright (c) 2016 ETCG
 *
 * License: MIT
 * Website: https://github.com/ETCG/pBar
 *
 * Look at the website for details on the license or
 * changes to the license itself.
 *
 *
 * The original class has been splitted into 'pBar.h' and
 * 'pBar.cpp' for easier including by Matteo Seclì
 * (<secli.matteo@gmail.com>)
 */

#include "pBar.h"

//Namespace:
using namespace std;
using namespace pBarNamespace;

void pBar::update(double newProgress) {
    currentProgress += newProgress;
    amountOfFiller = (int)((currentProgress / neededProgress)*(double)pBarLength);
}
void pBar::print() {
    currUpdateVal %= pBarUpdater.length();
    cout << "\r" //Bring cursor to start of line
         << firstPartOfpBar; //Print out first part of pBar
    for (int a = 0; a < amountOfFiller; a++) { //Print out current progress
        cout << pBarFiller;
    }
    cout << pBarUpdater[currUpdateVal];
    for (int b = 0; b < pBarLength - amountOfFiller; b++) { //Print out spaces
        cout << " ";
    }
    cout << lastPartOfpBar //Print out last part of progress bar
         << " (" << (int)(100*(currentProgress/neededProgress)) << "%)" //This just prints out the percent
         << flush;
    currUpdateVal += 1;
}
