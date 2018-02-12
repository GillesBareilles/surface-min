#ifndef FORMES_H_INCLUDED
#define FORMES_H_INCLUDED

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

class Forme
{
protected:
  int ordre_;
  std::vector<int> liste_noeuds_;

public:
  Forme() : ordre_(-1), liste_noeuds_(std::vector<int>()) {};
  Forme& operator= (const Forme&);
  void print() const;
  int operator()(int) const;
  int operator[](int) const;
  int ordre() const;
  int nb_noeuds() const;
};


class Point : public Forme
{
public:
  Point(int);
};


class Segment : public Forme
{
public:
  Segment(int i, int j);
  Segment(int i, int j, int k);
};


class Triangle : public  Forme
{
public:
  Triangle(int, int, int);
  Triangle(int, int, int, int, int, int);
};

#endif
