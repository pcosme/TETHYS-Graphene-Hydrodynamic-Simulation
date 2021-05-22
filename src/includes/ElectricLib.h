/************************************************************************************************\
* Copyright (c) 2020 Pedro Cosme and João Santos                                                 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

#ifndef ELECTRICLIB_H
#define ELECTRICLIB_H

#include <H5Cpp.h>
#include "TethysBaseLib.h"
#include "TethysMathLib.h"
#include "Fluid1DLib.h"
#include "Fluid2DLib.h"
#include "GrapheneFluid2DLib.h"
#include "GrapheneFluid1DLib.h"

/*!
 * @brief Class to obtain the macroscopic electrical quantities
 *
 * The main simulation computes the time evolution of the mean field quantities of number density and velocity.
 * From these the appropriate macroscopic quantities need to be derived in order to benchmark the numerical results with literature or experiments.
 * This class provides the methods to compute such quantities from the data stored on the GrapheneFluid1D or GrapheneFluid2D classes.
 *
 **/
class ElectroAnalysis{
private:
	std::ofstream data_electro;
	vector<float> TmpArr;
	vector<float> NetQ;
	vector<float> DipX;
	vector<float> DipY;
	vector<float> DipVarX;
	vector<float> DipVarY;
	vector<float> DipVarVarX;
	vector<float> DipVarVarY;
	vector<float> AvgCurDS;
	vector<float> VoltDS;
	vector<float> CurS;
	vector<float> CurD;
	vector<float> AvgCurHall;
	vector<float> PowOhm;
	vector<float> EngCap;
	vector<float> PowCap;

public:
	void CloseElectroFile();
	void CreateElectroFile(const GrapheneFluid1D& graphene);
	void CreateElectroFile(const string &infix_string);
	void WriteElectroFile(float t,const GrapheneFluid1D& graphene);
	void WriteElectroFile();

	/*!
	 * @brief 1D Integrator for the total electric charge.
	 *
	 *
	 * Computes the integral @f$ \int_0^1 n \,dx @f$
	 * */
	static float NetCharge(const GrapheneFluid1D& graphene);

	/*!
	 * @brief 2D Integrator for the total electric charge.
	 *
	 * The total charge is given by
	 * @f[ Q = e\int_0^W\int_0^L n\,dxdy = en_0L^2 \int_0^w\int_0^1 n^*\,dx^*dy^* @f] thus the charge prefactor is @f$ Q_0=en_0L^2 @f$. This method Computes the integral @f$ \int_0^{w}\int_0^1 n \,dxdy @f$
	 * */
	static float NetCharge(const GrapheneFluid2D& graphene);

	/*!
	 * @brief 1D Integrator for the dissipated Ohm power.
	 *
	 * This method computes the integral @f$ \int_0^1 jv \,dx @f$
	 * */
	static float OhmPower(const GrapheneFluid1D& graphene);

	/*!
	 * @brief 2D Integrator for the dissipated Ohm power.
	 *
	 * The dissipated Ohm power can be estimated with @f$ P_\Omega = J\cdot E @f$  nothing that @f$ J = \sigma E @f$  and @f$ \sigma=ne^2\tau/m^\star \doteq \sigma_0\sqrt{n^*}  @f$. Therefore,
	 *
	 * @f[ P_\Omega = \int_0^W\int_0^L \frac{j_x^2+j_y^2}{\sigma}\,dxdy = \frac{e^2v_0^2n_0^2L^2}{\sigma_0} \int_0^w\int_0^1 \frac{{j_x^*}^2+{j_y^*}^2}{\sqrt{n^*}}\,dx^*dy^*  @f]
	 *
	 * This method computes the integral @f$ \int_0^{w}\int_0^1 |j|^2/\sqrt{n} \,dxdy @f$
	 * */
	static float OhmPower(const GrapheneFluid2D& graphene);

	/*!
	* @brief 1D Integrator for the average drain-to-source current
	*
	* This method computes the integral @f$ \int_0^1 j \,dx @f$
	* */
	static float AverageCurrent(const GrapheneFluid1D& graphene);

	/*!
	 * @brief Integrator for the average Hall current
	 *
	 * The average Hall current is given by
     * @f[ I_{Hall}= en_0v_0L \frac{L}{W}\int_0^{w}\int_0 j_y^*\,dx^*dy^* @f]
	 * and the current prefactor is given by @f$ en_0v_0L@f$
	 *
	 * This method computes the integral @f$ \int_0^{w}\int_0^1 j_y \,dxdy @f$
	 * */
	static float AverageHallCurrent(const GrapheneFluid2D& graphene);

	/*!
	* @brief 2D Integrator for the average drain-to-source current
	*
	* The average  drain-to-source current @f$ I_{DS}=\frac{1}{L}\int_0^L I(x)\,dx @f$ with  @f$ I(x) = \int_0^W ej_x \,dy   @f$. Thus
	 *
	 * @f[ I_{DS}= en_0v_0L \int_0^{w}\int_0 j_x^*\,dx^*dy^* @f]
	* and the current prefactor is given by @f$ en_0v_0L@f$
	* This method computes the integral @f$ \int_0^{w}\int_0^1 j_x \,dxdy @f$
	* */
	static float AverageDirectCurrent(const GrapheneFluid2D& graphene);

	/*!
	* @brief Integrator for the local drain current
	*
	* This method computes the integral @f$ \int_0^{w} j_x(1,y) \,dy @f$
	* */
	static float DrainCurrent(const GrapheneFluid2D& graphene);

	/*!
	* @brief Integrator for the local source current
	*
	* This method computes the integral @f$ \int_0^{w} j_x(0,y) \,dy @f$
	* */
	static float SourceCurrent(const GrapheneFluid2D& graphene);

	/*!
	* @brief Integrator for the drain-to-source voltage
	*
	* This method computes the integral @f$ \int_0^{w}\int_0^1 v_x(1,y) \,dxdy @f$
	* */
	static float DrainToSourceVoltage(const GrapheneFluid2D& graphene);

	/*!
	* @brief 1D Integrator for the electric dipole moment
	*
	* This method computes the integral @f$ \int_0^1 (x-1/2)n \,dx @f$
	* */
	static float ElectricDipole(const GrapheneFluid1D& graphene);

	/*!
	* @brief 2D Integrator for the electric dipole moment x component
	*
	* This method computes the integral @f$ \int_0^{W/L}\int_0^1 (x-1/2)n \,dxdy @f$
	* */
	static float ElectricDipoleX(const GrapheneFluid2D& graphene);

	/*!
	* @brief  2D Integrator for the electric dipole moment y component
	*
	* This method computes the integral @f$ \int_0^{w}\int_0^1 (y-w/2)n \,dxdy @f$
	* */
	static float ElectricDipoleY(const GrapheneFluid2D& graphene);

	/*!
	* @brief 1D Integrator for the electric dipole moment derivative
	*
	* This method computes the integral @f$ \int_0^1 j \,dx @f$
	* */
	static float ElectricDipoleVariation(const GrapheneFluid1D& graphene);

	/*!
	* @brief 2D Integrator for the electric dipole moment derivative x component
	*
	* This method computes the integral @f$ \int_0^{w}\int_0^1 j_x \,dxdy @f$
	* */
	static float ElectricDipoleVariationX(const GrapheneFluid2D& graphene);

	/*!
	* @brief 2D Integrator for the electric dipole moment derivative y component
	*
	* This method computes the integral @f$ \int_0^{w}\int_0^1 j_y \,dxdy @f$
	* */
	static float ElectricDipoleVariationY(const GrapheneFluid2D& graphene);


	/*!
	 * @brief Computes the electric quantities that can be obtained by direct integration.
	 *
	 * This method computes directly the following quantities
	 * - Net electric charge
	 * - Electric dipole moment
	 * - Local current at drain
	 * - Local current at source
	 * - Average drain-to-source current
	 * - Average Hall current
	 * - Drain-to-source voltage
	 * - Dissipated Ohm Power
	 *
	 * @param t time of the snapshot
	 * @param graphene GrapheneFluid2D instance from where the quantities are computed
	 * */
	void ComputeElectroBase(float t, const GrapheneFluid2D &graphene);

	/*!
	 * @brief Computes the electric quantities that can not be obtained by direct integration.
	 *
	 * This method computes directly the following quantities
	 * - Stored power at the gate/graphene capacitor
	 * - First time derivative of the electric dipole momement
	 * - Second time derivative of the electric dipole momement
	 *
	 * @param t time of the snapshot
	 * @param graphene GrapheneFluid2D instance from where the quantities are computed
	 * */
	void ComputeElectroDerived();

	static void BannerDisplay(const GrapheneFluid2D &graphene); ///< @brief Displays the electrical simulations banner
};

#endif
