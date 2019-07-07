
#ifndef HELPERS_H
#define HELPERS_H

#include <math.h>
#include <string>
#include <vector>
#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>

using namespace std;

constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

//====
// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
	auto found_null = s.find("null");
	auto b1 = s.find_first_of("[");
	auto b2 = s.find_first_of("}");
	if (found_null != string::npos) {
		return "";
	}
	else if (b1 != string::npos && b2 != string::npos) {
		return s.substr(b1, b2 - b1 + 2);
	}
	return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for (int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x, y, map_x, map_y);
		if (dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x, y, maps_x, maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y - y), (map_x - x));

	double angle = abs(theta - heading);

	if (angle > pi() / 4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x, y, theta, maps_x, maps_y);

	int prev_wp;
	prev_wp = next_wp - 1;
	if (next_wp == 0)
	{
		prev_wp = maps_x.size() - 1;
	}

	double n_x = maps_x[next_wp] - maps_x[prev_wp];
	double n_y = maps_y[next_wp] - maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x + x_y*n_y) / (n_x*n_x + n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x, x_y, proj_x, proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000 - maps_x[prev_wp];
	double center_y = 2000 - maps_y[prev_wp];
	double centerToPos = distance(center_x, center_y, x_x, x_y);
	double centerToRef = distance(center_x, center_y, proj_x, proj_y);

	if (centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for (int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i], maps_y[i], maps_x[i + 1], maps_y[i + 1]);
	}

	frenet_s += distance(0, 0, proj_x, proj_y);

	return{ frenet_s,frenet_d };

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while (s > maps_s[prev_wp + 1] && (prev_wp < (int)(maps_s.size() - 1)))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp + 1) % maps_x.size();

	double heading = atan2((maps_y[wp2] - maps_y[prev_wp]), (maps_x[wp2] - maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s - maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp] + seg_s*cos(heading);
	double seg_y = maps_y[prev_wp] + seg_s*sin(heading);

	double perp_heading = heading - pi() / 2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return{ x,y };

}

namespace
{

	namespace hp
	{

		// band matrix solver
		class band_matrix
		{
		private:
			std::vector< std::vector<double> > m_upper;  // upper band
			std::vector< std::vector<double> > m_lower;  // lower band
		public:
			band_matrix() {};                             // constructor
			band_matrix(int dim, int n_u, int n_l);       // constructor
			~band_matrix() {};                            // destructor
			void resize(int dim, int n_u, int n_l);      // init with dim,n_u,n_l
			int dim() const;                             // matrix dimension
			int num_upper() const
			{
				return m_upper.size() - 1;
			}
			int num_lower() const
			{
				return m_lower.size() - 1;
			}
			// access operator
			double & operator () (int i, int j);            // write
			double   operator () (int i, int j) const;      // read
															// we can store an additional diogonal (in m_lower)
			double& saved_diag(int i);
			double  saved_diag(int i) const;
			void lu_decompose();
			std::vector<double> r_solve(const std::vector<double>& b) const;
			std::vector<double> l_solve(const std::vector<double>& b) const;
			std::vector<double> lu_solve(const std::vector<double>& b,
				bool is_lu_decomposed = false);

		};

		band_matrix::band_matrix(int dim, int n_u, int n_l)
		{
			resize(dim, n_u, n_l);
		}
		void band_matrix::resize(int dim, int n_u, int n_l)
		{
			assert(dim>0);
			assert(n_u >= 0);
			assert(n_l >= 0);
			m_upper.resize(n_u + 1);
			m_lower.resize(n_l + 1);
			for (size_t i = 0; i<m_upper.size(); i++) {
				m_upper[i].resize(dim);
			}
			for (size_t i = 0; i<m_lower.size(); i++) {
				m_lower[i].resize(dim);
			}
		}
		int band_matrix::dim() const
		{
			if (m_upper.size()>0) {
				return m_upper[0].size();
			}
			else {
				return 0;
			}
		}


		// defines the new operator (), so that we can access the elements
		// by A(i,j), index going from i=0,...,dim()-1
		double & band_matrix::operator () (int i, int j)
		{
			int k = j - i;       // what band is the entry
			assert((i >= 0) && (i<dim()) && (j >= 0) && (j<dim()));
			assert((-num_lower() <= k) && (k <= num_upper()));
			// k=0 -> diogonal, k<0 lower left part, k>0 upper right part
			if (k >= 0)   return m_upper[k][i];
			else	    return m_lower[-k][i];
		}
		double band_matrix::operator () (int i, int j) const
		{
			int k = j - i;       // what band is the entry
			assert((i >= 0) && (i<dim()) && (j >= 0) && (j<dim()));
			assert((-num_lower() <= k) && (k <= num_upper()));
			// k=0 -> diogonal, k<0 lower left part, k>0 upper right part
			if (k >= 0)   return m_upper[k][i];
			else	    return m_lower[-k][i];
		}
		// second diag (used in LU decomposition), saved in m_lower
		double band_matrix::saved_diag(int i) const
		{
			assert((i >= 0) && (i<dim()));
			return m_lower[0][i];
		}
		double & band_matrix::saved_diag(int i)
		{
			assert((i >= 0) && (i<dim()));
			return m_lower[0][i];
		}

		// LR-Decomposition of a band matrix
		void band_matrix::lu_decompose()
		{
			int  i_max, j_max;
			int  j_min;
			double x;

			// preconditioning
			// normalize column i so that a_ii=1
			for (int i = 0; i<this->dim(); i++) {
				assert(this->operator()(i, i) != 0.0);
				this->saved_diag(i) = 1.0 / this->operator()(i, i);
				j_min = std::max(0, i - this->num_lower());
				j_max = std::min(this->dim() - 1, i + this->num_upper());
				for (int j = j_min; j <= j_max; j++) {
					this->operator()(i, j) *= this->saved_diag(i);
				}
				this->operator()(i, i) = 1.0;          // prevents rounding errors
			}

			// Gauss LR-Decomposition
			for (int k = 0; k<this->dim(); k++) {
				i_max = std::min(this->dim() - 1, k + this->num_lower());  // num_lower not a mistake!
				for (int i = k + 1; i <= i_max; i++) {
					assert(this->operator()(k, k) != 0.0);
					x = -this->operator()(i, k) / this->operator()(k, k);
					this->operator()(i, k) = -x;                         // assembly part of L
					j_max = std::min(this->dim() - 1, k + this->num_upper());
					for (int j = k + 1; j <= j_max; j++) {
						// assembly part of R
						this->operator()(i, j) = this->operator()(i, j) + x*this->operator()(k, j);
					}
				}
			}
		}
		// solves Ly=b
		std::vector<double> band_matrix::l_solve(const std::vector<double>& b) const
		{
			assert(this->dim() == (int)b.size());
			std::vector<double> x(this->dim());
			int j_start;
			double sum;
			for (int i = 0; i<this->dim(); i++) {
				sum = 0;
				j_start = std::max(0, i - this->num_lower());
				for (int j = j_start; j<i; j++) sum += this->operator()(i, j)*x[j];
				x[i] = (b[i] * this->saved_diag(i)) - sum;
			}
			return x;
		}
		// solves Rx=y
		std::vector<double> band_matrix::r_solve(const std::vector<double>& b) const
		{
			assert(this->dim() == (int)b.size());
			std::vector<double> x(this->dim());
			int j_stop;
			double sum;
			for (int i = this->dim() - 1; i >= 0; i--) {
				sum = 0;
				j_stop = std::min(this->dim() - 1, i + this->num_upper());
				for (int j = i + 1; j <= j_stop; j++) sum += this->operator()(i, j)*x[j];
				x[i] = (b[i] - sum) / this->operator()(i, i);
			}
			return x;
		}

		std::vector<double> band_matrix::lu_solve(const std::vector<double>& b,
			bool is_lu_decomposed)
		{
			assert(this->dim() == (int)b.size());
			std::vector<double>  x, y;
			if (is_lu_decomposed == false) {
				this->lu_decompose();
			}
			y = this->l_solve(b);
			x = this->r_solve(y);
			return x;
		}

		class spline
		{
		public:
			enum bd_type {
				first_deriv = 1,
				second_deriv = 2
			};

		private:
			std::vector<double> m_x, m_y;            // x,y coordinates of points
													 // interpolation parameters
													 // f(x) = a*(x-x_i)^3 + b*(x-x_i)^2 + c*(x-x_i) + y_i
			std::vector<double> m_a, m_b, m_c;        // spline coefficients
			double  m_b0, m_c0;                     // for left extrapol
			bd_type m_left, m_right;
			double  m_left_value, m_right_value;
			bool    m_force_linear_extrapolation;

		public:
			// set default boundary condition to be zero curvature at both ends
			spline() : m_left(second_deriv), m_right(second_deriv),
				m_left_value(0.0), m_right_value(0.0),
				m_force_linear_extrapolation(false)
			{
				;
			}

			// optional, but if called it has to come be before set_points()
			void set_boundary(bd_type left, double left_value,
				bd_type right, double right_value,
				bool force_linear_extrapolation = false);
			void set_points(const std::vector<double>& x,
				const std::vector<double>& y, bool cubic_spline = true);
			double operator() (double x) const;
		};

		void spline::set_boundary(spline::bd_type left, double left_value,
			spline::bd_type right, double right_value,
			bool force_linear_extrapolation)
		{
			assert(m_x.size() == 0);          // set_points() must not have happened yet
			m_left = left;
			m_right = right;
			m_left_value = left_value;
			m_right_value = right_value;
			m_force_linear_extrapolation = force_linear_extrapolation;
		}


		void spline::set_points(const std::vector<double>& x,
			const std::vector<double>& y, bool cubic_spline)
		{
			assert(x.size() == y.size());
			assert(x.size()>2);
			m_x = x;
			m_y = y;
			int   n = x.size();
			// TODO: maybe sort x and y, rather than returning an error
			for (int i = 0; i<n - 1; i++) {
				assert(m_x[i]<m_x[i + 1]);
			}

			if (cubic_spline == true) { 
				band_matrix A(n, 1, 1);
				std::vector<double>  rhs(n);
				for (int i = 1; i<n - 1; i++) {
					A(i, i - 1) = 1.0 / 3.0*(x[i] - x[i - 1]);
					A(i, i) = 2.0 / 3.0*(x[i + 1] - x[i - 1]);
					A(i, i + 1) = 1.0 / 3.0*(x[i + 1] - x[i]);
					rhs[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
				}
				// boundary conditions
				if (m_left == spline::second_deriv) {
					// 2*b[0] = f''
					A(0, 0) = 2.0;
					A(0, 1) = 0.0;
					rhs[0] = m_left_value;
				}
				else if (m_left == spline::first_deriv) {
					// c[0] = f', needs to be re-expressed in terms of b:
					// (2b[0]+b[1])(x[1]-x[0]) = 3 ((y[1]-y[0])/(x[1]-x[0]) - f')
					A(0, 0) = 2.0*(x[1] - x[0]);
					A(0, 1) = 1.0*(x[1] - x[0]);
					rhs[0] = 3.0*((y[1] - y[0]) / (x[1] - x[0]) - m_left_value);
				}
				else {
					assert(false);
				}
				if (m_right == spline::second_deriv) {
					// 2*b[n-1] = f''
					A(n - 1, n - 1) = 2.0;
					A(n - 1, n - 2) = 0.0;
					rhs[n - 1] = m_right_value;
				}
				else if (m_right == spline::first_deriv) {
					// c[n-1] = f', needs to be re-expressed in terms of b:
					// (b[n-2]+2b[n-1])(x[n-1]-x[n-2])
					// = 3 (f' - (y[n-1]-y[n-2])/(x[n-1]-x[n-2]))
					A(n - 1, n - 1) = 2.0*(x[n - 1] - x[n - 2]);
					A(n - 1, n - 2) = 1.0*(x[n - 1] - x[n - 2]);
					rhs[n - 1] = 3.0*(m_right_value - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
				}
				else {
					assert(false);
				}

				// solve the equation system to obtain the parameters b[]
				m_b = A.lu_solve(rhs);

				// calculate parameters a[] and c[] based on b[]
				m_a.resize(n);
				m_c.resize(n);
				for (int i = 0; i<n - 1; i++) {
					m_a[i] = 1.0 / 3.0*(m_b[i + 1] - m_b[i]) / (x[i + 1] - x[i]);
					m_c[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
						- 1.0 / 3.0*(2.0*m_b[i] + m_b[i + 1])*(x[i + 1] - x[i]);
				}
			}
			else { // linear interpolation
				m_a.resize(n);
				m_b.resize(n);
				m_c.resize(n);
				for (int i = 0; i<n - 1; i++) {
					m_a[i] = 0.0;
					m_b[i] = 0.0;
					m_c[i] = (m_y[i + 1] - m_y[i]) / (m_x[i + 1] - m_x[i]);
				}
			}

			// for left extrapolation coefficients
			m_b0 = (m_force_linear_extrapolation == false) ? m_b[0] : 0.0;
			m_c0 = m_c[0];

			// for the right extrapolation coefficients
			// f_{n-1}(x) = b*(x-x_{n-1})^2 + c*(x-x_{n-1}) + y_{n-1}
			double h = x[n - 1] - x[n - 2];
			// m_b[n-1] is determined by the boundary condition
			m_a[n - 1] = 0.0;
			m_c[n - 1] = 3.0*m_a[n - 2] * h*h + 2.0*m_b[n - 2] * h + m_c[n - 2];   // = f'_{n-2}(x_{n-1})
			if (m_force_linear_extrapolation == true)
				m_b[n - 1] = 0.0;
		}

		double spline::operator() (double x) const
		{
			size_t n = m_x.size();
			// find the closest point m_x[idx] < x, idx=0 even if x<m_x[0]
			std::vector<double>::const_iterator it;
			it = std::lower_bound(m_x.begin(), m_x.end(), x);
			int idx = std::max(int(it - m_x.begin()) - 1, 0);

			double h = x - m_x[idx];
			double interpol;
			if (x<m_x[0]) {
				// extrapolation to the left
				interpol = (m_b0*h + m_c0)*h + m_y[0];
			}
			else if (x>m_x[n - 1]) {
				// extrapolation to the right
				interpol = (m_b[n - 1] * h + m_c[n - 1])*h + m_y[n - 1];
			}
			else {
				// interpolation
				interpol = ((m_a[idx] * h + m_b[idx])*h + m_c[idx])*h + m_y[idx];
			}
			return interpol;
		}


	} // namespace hp


} // namespace

#endif 