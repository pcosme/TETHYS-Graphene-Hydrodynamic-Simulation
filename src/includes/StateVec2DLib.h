/************************************************************************************************\
* 2022 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* 																 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


#ifndef STATEVECLIB1D_H
#define STATEVECLIB1D_H

#include <iostream>

class StateVec2D{
private:
	float density;
	float velocity_x;
	float velocity_y;
	float temperature;
public:
	StateVec2D()=default; //default constructor
	StateVec2D(const StateVec2D&); //copy constructor
	StateVec2D(float den, float velx, float vely,float temp);
	~StateVec2D()=default;

	float& vx();
	float& vy();
	float& n();
	float& tmp();

	StateVec2D& operator=(const StateVec2D&);
	StateVec2D operator + (StateVec2D const &obj) const ;
	StateVec2D operator - (StateVec2D const &obj) const ;
	StateVec2D operator * (StateVec2D const &obj) const ;
	StateVec2D operator / (StateVec2D const &obj) const ;
	friend StateVec2D operator*(const StateVec2D &obj, float value);
	friend StateVec2D operator*(float value, const StateVec2D &obj);
	friend StateVec2D operator/(const StateVec2D &obj, float value);
	friend std::ostream& operator<<(std::ostream& outstream, const StateVec2D &obj);
};




#endif //STATEVECLIB1D_H


