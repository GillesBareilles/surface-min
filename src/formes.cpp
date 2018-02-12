#include "formes.hpp"

//---------------------- Forme
Forme& Forme::operator= (const Forme& F) {
  if (liste_noeuds_==this->liste_noeuds_) return *this;
  ordre_ = F.ordre_;
  // nb_noeuds_ = F.nb_noeuds_;
  std::vector<int> liste_noeuds_(F.liste_noeuds_.size());
  for (unsigned int i=0; i<F.liste_noeuds_.size(); i++) liste_noeuds_[i] = F.liste_noeuds_[i];
  return *this;
}

int Forme::operator()(int i) const {
  return liste_noeuds_.at(i-1);
}

int Forme::operator[](int i) const {
  return liste_noeuds_.at(i);
}

void Forme::print() const {
  for (int i:liste_noeuds_) std::cout << i << " ";
  std::cout << std::endl;
}

int Forme::ordre() const {
  return ordre_;
}

int Forme::nb_noeuds() const {
  return liste_noeuds_.size();
}

//---------------------- Point
Point::Point(int i) : Forme()
{
    ordre_ = 0;
    liste_noeuds_.push_back(i);
}

//---------------------- Segment
Segment::Segment(int i, int j) : Forme() {
    ordre_ = 1;
    liste_noeuds_.push_back(i);
    liste_noeuds_.push_back(j);
}

Segment::Segment(int i, int j, int k) : Forme() {
    ordre_ = 2;
    liste_noeuds_.push_back(i);
    liste_noeuds_.push_back(j);
    liste_noeuds_.push_back(k);
}


//---------------------- Triangle
Triangle::Triangle(int i, int j, int k) : Forme() {
    ordre_ = 1;
    liste_noeuds_.push_back(i);
    liste_noeuds_.push_back(j);
    liste_noeuds_.push_back(k);
}

Triangle::Triangle(int i, int j, int k, int l, int m, int n) : Forme() {
    ordre_ = 2;
    liste_noeuds_.push_back(i);
    liste_noeuds_.push_back(j);
    liste_noeuds_.push_back(k);
    liste_noeuds_.push_back(l);
    liste_noeuds_.push_back(m);
    liste_noeuds_.push_back(n);
}
