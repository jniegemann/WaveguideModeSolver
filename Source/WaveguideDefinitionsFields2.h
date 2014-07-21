#ifndef _WAVEGUIDEFIELDS_H
#define _WAVEGUIDEFIELDS_H

#include "hermes2d.h"

template<typename T>
struct Material
{
  std::string name;

  T epsXX;
  T epsYY;
  T epsZZ;

  T invMuXX;
  T invMuYY;
  T invMuZZ;

  Material(const std::string& in_name, const T& in_epsXX, const T& in_epsYY, const T& in_epsZZ, const T& in_invMuXX, const T& in_invMuYY, const T& in_invMuZZ)
  : name(in_name), epsXX(in_epsXX), epsYY(in_epsYY), epsZZ(in_epsZZ), invMuXX(in_invMuXX), invMuYY(in_invMuYY), invMuZZ(in_invMuZZ) {}
};


/**
 * A custom bilinear form A
 */
template<typename T>
class CustomMatrixA : public Hermes::Hermes2D::MatrixFormVol<T>   
{
private:
  const double k0;
  const T epsXX;
  const T epsYY;
  const T invMuZZ;

  using Hermes::Hermes2D::Form<T>::areas;
  using Hermes::Hermes2D::MatrixForm<T>::i;
  using Hermes::Hermes2D::MatrixForm<T>::j;
public:
  
  CustomMatrixA(int i, int j, std::string in_area, double in_k0, T in_epsXX, T in_epsYY, T in_invMuZZ)
  : Hermes::Hermes2D::MatrixFormVol<T>(i, j), k0(in_k0), epsXX(in_epsXX), epsYY(in_epsYY), invMuZZ(in_invMuZZ)
  {
    this->set_area(in_area);
    this->setSymFlag(Hermes::Hermes2D::HERMES_NONSYM);
  }

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<Real> *u, Hermes::Hermes2D::Func<Real> *v, Hermes::Hermes2D::Geom<Real> *e, Hermes::Hermes2D::Func<Scalar>  **ext) const
  {
    Scalar result(0);
    for (int i = 0; i < n; i++)
    {
      result += wt[i] * ((u->curl[i] * invMuZZ * v->curl[i]) - k0*k0*(u->val0[i] * epsXX * v->val0[i] + u->val1[i] * epsYY * v->val1[i]));
    }

    return result;
  };

  virtual T value(int n, double *wt, Hermes::Hermes2D::Func<T> *u_ext[], Hermes::Hermes2D::Func<double> *u, Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<T>  **ext) const
  {
    return matrix_form<double, T>(n, wt, u_ext, u, v, e, ext);
  }

  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u, Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord>  **ext) const
  {
    return matrix_form<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, u, v, e, ext);
  };

  Hermes::Hermes2D::MatrixFormVol<T>* clone() const
  {
    return new CustomMatrixA(i, j, areas[0], k0, epsXX, epsYY, invMuZZ);
  };

};


/*
 * A custom bilinear form B
 */
template<typename T>
class CustomMatrixB : public Hermes::Hermes2D::MatrixFormVol<T>   
{
private:
  const T invMuXX;
  const T invMuYY;
  using Hermes::Hermes2D::Form<T>::areas;
  using Hermes::Hermes2D::MatrixForm<T>::i;
  using Hermes::Hermes2D::MatrixForm<T>::j;
public:
  
  CustomMatrixB(int i, int j, std::string in_area, T in_invMuXX, T in_invMuYY)
  : Hermes::Hermes2D::MatrixFormVol<T>(i, j), invMuXX(in_invMuXX), invMuYY(in_invMuYY)
  {
    this->set_area(in_area);
    this->setSymFlag(Hermes::Hermes2D::HERMES_NONSYM);
  }

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<Real> *u, Hermes::Hermes2D::Func<Real> *v, Hermes::Hermes2D::Geom<Real> *e, Hermes::Hermes2D::Func<Scalar>  **ext) const
  {
    Scalar result(0);
    for (int i = 0; i < n; i++)
    {
      result += wt[i] * (u->val0[i] * invMuYY * v->val0[i] + u->val1[i] * invMuXX * v->val1[i]);
    }

    return result;
  };

  virtual T value(int n, double *wt, Hermes::Hermes2D::Func<T> *u_ext[], Hermes::Hermes2D::Func<double> *u, Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<T>  **ext) const
  {
    return matrix_form<double, T>(n, wt, u_ext, u, v, e, ext);
  }

  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u, Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord>  **ext) const
  {
    return matrix_form<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, u, v, e, ext);
  };

  Hermes::Hermes2D::MatrixFormVol<T>* clone() const
  {
    return new CustomMatrixB(i, j, areas[0], invMuXX, invMuYY);
  };

};

/**
 * A custom bilinear form for the mixed term grad(phi) dot v, i.e. the transpose of CustomMatrixFormMixedZT
 */
template<typename T>
class CustomMatrixC : public Hermes::Hermes2D::MatrixFormVol<T>   
{
private:
  const T invMuXX;
  const T invMuYY;

  using Hermes::Hermes2D::Form<T>::areas;
  using Hermes::Hermes2D::MatrixForm<T>::i;
  using Hermes::Hermes2D::MatrixForm<T>::j;
public:
  
  CustomMatrixC(int i, int j, std::string in_area, const T& in_invMuXX, const T& in_invMuYY)
  : Hermes::Hermes2D::MatrixFormVol<T>(i, j), invMuXX(in_invMuXX), invMuYY(in_invMuYY)
  {
    this->set_area(in_area);
    this->setSymFlag(Hermes::Hermes2D::HERMES_NONSYM);
  }
 
  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<Real> *u, Hermes::Hermes2D::Func<Real> *v, Hermes::Hermes2D::Geom<Real> *e, Hermes::Hermes2D::Func<Scalar>  **ext) const
  {
    Scalar result = Scalar(0); 
    for (int i = 0; i < n; i++) result += wt[i] * (u->dx[i] * invMuYY * v->val0[i] + u->dy[i] * invMuXX * v->val1[i]);
    return result;
  };

  virtual T value(int n, double *wt, Hermes::Hermes2D::Func<T> *u_ext[], Hermes::Hermes2D::Func<double> *u, Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<T>  **ext) const  {
    return matrix_form<double, T>(n, wt, u_ext, u, v, e, ext);
  }

  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u, Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord>  **ext) const
  {
    return matrix_form<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, u, v, e, ext);
  };

  Hermes::Hermes2D::MatrixFormVol<T>* clone() const
  {
    return new CustomMatrixC(i, j, areas[0], invMuXX, invMuYY);
  };

private:

  T const_coeff;
  Hermes::Hermes2DFunction<T>* function_coeff;
  Hermes::Hermes2D::GeomType gt;
};

/**
 * A custom bilinear form D
 */
template<typename T>
class CustomMatrixD : public Hermes::Hermes2D::MatrixFormVol<T>   
{
private:
  const double k0;
  const T invMuXX;
  const T invMuYY;
  const T epsZZ;

  using Hermes::Hermes2D::Form<T>::areas;
  using Hermes::Hermes2D::MatrixForm<T>::i;
  using Hermes::Hermes2D::MatrixForm<T>::j;
public:
  
  CustomMatrixD(int i, int j, std::string in_area, double in_k0, T in_epsZZ, T in_invMuXX, T in_invMuYY)
  : Hermes::Hermes2D::MatrixFormVol<T>(i, j), k0(in_k0), epsZZ(in_epsZZ), invMuXX(in_invMuXX), invMuYY(in_invMuYY)
  {
    this->set_area(in_area);
    this->setSymFlag(Hermes::Hermes2D::HERMES_NONSYM);
  }

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<Real> *u, Hermes::Hermes2D::Func<Real> *v, Hermes::Hermes2D::Geom<Real> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
  {
    Scalar result(0);
    for (int i = 0; i < n; i++)
    {
      result += wt[i] * ((u->dx[i] * invMuYY * v->dx[i] + u->dy[i] * invMuXX * v->dy[i]) - k0*k0*(u->val[i] * epsZZ * v->val[i]) );
   }

    return result;
  };

  virtual T value(int n, double *wt, Hermes::Hermes2D::Func<T> *u_ext[], Hermes::Hermes2D::Func<double> *u, Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<T>  **ext) const
  {
    return matrix_form<double, T>(n, wt, u_ext, u, v, e, ext);
  }

  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u, Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord>  **ext) const
  {
    return matrix_form<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, u, v, e, ext);
  };

  Hermes::Hermes2D::MatrixFormVol<T>* clone() const
  {
    return new CustomMatrixD(i, j, areas[0], k0, epsZZ, invMuXX, invMuYY);
  };

};



/**
 * A custom bilinear form for the mixed term v dot grad(phi) where v is the transversal vector field and phi the scalar z-component
 */
template<typename T>
class CustomMatrixE : public Hermes::Hermes2D::MatrixFormVol<T>   
{
private:
  const T invMuXX;
  const T invMuYY;

  using Hermes::Hermes2D::Form<T>::areas;
  using Hermes::Hermes2D::MatrixForm<T>::i;
  using Hermes::Hermes2D::MatrixForm<T>::j;
public:
  
  CustomMatrixE(int i, int j, std::string in_area, const T& in_invMuXX, const T& in_invMuYY)
  : Hermes::Hermes2D::MatrixFormVol<T>(i, j), invMuXX(in_invMuXX), invMuYY(in_invMuYY)
  {
    this->set_area(in_area);
    this->setSymFlag(Hermes::Hermes2D::HERMES_NONSYM);
  }

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<Real> *u, Hermes::Hermes2D::Func<Real> *v, Hermes::Hermes2D::Geom<Real> *e, Hermes::Hermes2D::Func<Scalar>  **ext) const
  {
    Scalar result = Scalar(0); 
    for (int i = 0; i < n; i++) result += wt[i] * (u->val0[i] * invMuYY * v->dx[i] + u->val1[i] * invMuXX * v->dy[i]);
    return result;
  };

  virtual T value(int n, double *wt, Hermes::Hermes2D::Func<T> *u_ext[], Hermes::Hermes2D::Func<double> *u, Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<T>  **ext) const
  {
    return matrix_form<double, T>(n, wt, u_ext, u, v, e, ext);
  }

  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u, Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord>  **ext) const
  {
    return matrix_form<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, u, v, e, ext);
  };

  CustomMatrixE<T>* clone() const
  {
    return new CustomMatrixE(i, j, areas[0], invMuXX, invMuYY);
  };

};







// This is the left-hand side of the eigenvalue problem for the waveguide modes, following the notation of Davidson (eq. 10.150)
template<typename T>
class WeakFormWaveguideLeft : public Hermes::Hermes2D::WeakForm<T>
{
        const double m_k0;
        const std::vector<Material<T>>& m_materialList;

	using Hermes::Hermes2D::WeakForm<T>::add_matrix_form;
public:
	WeakFormWaveguideLeft(double k0, const std::vector<Material<T>>& materialList) : Hermes::Hermes2D::WeakForm<T>(2), m_k0(k0), m_materialList(materialList)
	{
	  // Iterate over all materials in our system
	  for( const auto& curMaterial : materialList )
	  {
//		std::cout << "Adding material: " << curMaterial.name << std::endl;
	    // The following represents the A matrix
	    add_matrix_form(new CustomMatrixA<T>(0, 0, curMaterial.name, k0, curMaterial.epsXX, curMaterial.epsYY, curMaterial.invMuZZ) );
	    add_matrix_form(new CustomMatrixC<T>(0, 1, curMaterial.name, curMaterial.invMuXX, curMaterial.invMuYY));
	    add_matrix_form(new CustomMatrixD<T>(1, 1, curMaterial.name, k0, curMaterial.epsZZ, curMaterial.invMuXX, curMaterial.invMuYY) );
	  }
	}

	Hermes::Hermes2D::WeakForm<T>* clone() const
	{
	  return new WeakFormWaveguideLeft<T>(m_k0, m_materialList );
	}
};

// This is the right-hand  side of the eigenvalue problem for the waveguide modes, following the notation of Davidson (eq. 10.150)
template<typename T>
class WeakFormWaveguideRight : public Hermes::Hermes2D::WeakForm<T>
{
        const double m_k0;
        const std::vector<Material<T>>& m_materialList;

	using Hermes::Hermes2D::WeakForm<T>::add_matrix_form;
public:
	WeakFormWaveguideRight(double k0, const std::vector<Material<T>>& materialList) : Hermes::Hermes2D::WeakForm<T>(2), m_k0(k0), m_materialList(materialList)
	{
	  // Iterate over all materials in our system
	  for( const auto& curMaterial : materialList )
	  {
//		std::cout << "Adding material: " << curMaterial.name << std::endl;
	    // The following represents the B matrix
	    add_matrix_form(new CustomMatrixB<T>(0, 0, curMaterial.name, curMaterial.invMuXX, curMaterial.invMuYY) );
	    add_matrix_form(new CustomMatrixE<T>(1, 0, curMaterial.name, curMaterial.invMuXX, curMaterial.invMuYY));
	  }
	}

	Hermes::Hermes2D::WeakForm<T>* clone() const
	{
	  return new WeakFormWaveguideRight<T>(m_k0, m_materialList );
	}
};

#endif

