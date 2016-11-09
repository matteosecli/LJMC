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

#ifndef PBAR_H
#define PBAR_H

#include <iostream> //For I/O
#include <string> //For strings

// Define a new pBar namespace
namespace pBarNamespace
{

//Class:
class pBar {
public:
    void update(double newProgress);
    void print();
    std::string firstPartOfpBar = "[", //Change these at will (that is why I made them public)
    lastPartOfpBar = "]",
    pBarFiller = "|",
    pBarUpdater = "/-\\|";
private:
    int amountOfFiller,
    pBarLength = 50, //I would recommend NOT changing this
    currUpdateVal = 0; //Do not change
    double currentProgress = 0, //Do not change
    neededProgress = 100; //I would recommend NOT changing this
};

}

#endif // PBAR_H
