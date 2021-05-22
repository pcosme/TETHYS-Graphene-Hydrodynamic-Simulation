
#ifndef GRID2DLIB_H
#define GRID2DLIB_H

#include <cstdlib>

class GridPoint{
private:
	int size_x;
	int size_y;
	int stride;
public:
	GridPoint(int pos,int xpoints, int ypoints,bool mid);
	~GridPoint() = default;
	bool IsMidGrid = false;
	int NW; //primeiros vizinhos portanto pertencentes a outro tipo de grelha
	int NE;
	int SW;
	int SE;
	int NW2; //terceiros vizinhos mesmo tipo de grelha
	int NE2;
	int SW2;
	int SE2;
	int N; //segundos vizinhos mesmo tipo de grelha
	int S;
	int E;
	int W;
	int C; //ponto de referencia da grelha
	void FirstNeighbours();
	void SecondNeighbours();
	void ThirdNeighbours();
};


#endif //GRID2DLIB_H
