#include "CovMatHandler.h"

/**
 * Constructor
 */
CovMatHandler::CovMatHandler()
{}

/**
 * Resizing the tensor @covmat to the number of bins
 * @num_points: number of bins
 */
void CovMatHandler::set_num_points(size_t num_points)
{
	//resizing
	covmat.resize(num_points);
}

/**
 * This function reads the covariance matrices from file and makes them accessible for later usage
 * @filename: the base of the names of the covariance matrices; the rest should only be a number
 * @index: bin number
 * @ifile: object for reading a file
 * @line: string that temporary stores a line from the file
 * @tmp: vector that temporary stores a line from the file as doubles
 * @tmp_mat: matrix that temporary stores the @tmp and later adds them to @covmat
 * @pos_begin, @pos_end: helper for splitting up @line to extract the numbers
 */
void CovMatHandler::read_covmat(string filename, size_t index)
{
	ifstream ifile;
	string line;
	vector<double> tmp;
	vector<vector<double>> tmp_mat;
	size_t pos_begin, pos_end;
	
	//open the file
	ifile.open(filename);
	if(ifile.is_open())
	{
		//walk over the file line wise
		while(getline(ifile, line))
		{
			pos_begin = 0;
			//splitting up the first number in the line
			//this needs to be in an extra part beacause of the entry at position 0 of the line
			if(line.find(" ", pos_begin + 1) < line.size())
			{
				pos_end = line.find(" ", pos_begin + 1);
				tmp.push_back(atof(line.substr(0, pos_end - pos_begin).c_str()));
				pos_begin = pos_end;
			}
			
			//splitting up the rest of the line
			while(line.find(" ", pos_begin + 1) < line.size())
			{
				pos_end = line.find(" ", pos_begin + 1);
				tmp.push_back(atof(line.substr(pos_begin + 1, pos_end - pos_begin - 1).c_str()));
				pos_begin = pos_end;
			}
			
			//adding the line to the matrix and emptying it
			tmp_mat.push_back(tmp);
			tmp.clear();
			
		}	
		
		//adding the matrix to @covmat
		covmat[index] = tmp_mat;
		
		is_symm(index);
		
		//close the file
		ifile.close();
		
		//construct a triangle matrix that can be computed faster
		make_triangle(index);
		
	}
	else
		cout << "error: could not read " << filename << endl;
	
}

/**
 * This function sets up the template matrix. This matrix contains the power of every parameter in every entry.
 * @power: list of all powers
 */
void CovMatHandler::set_templatemat(vector<vector<int>>& power)
{
	//start the first resizing to get the number of rows
	templatemat.resize(power.size());
	
	for(size_t row = 0; row < power.size(); row++)
	{
		//resizing for the number of coloumns
		templatemat[row].resize(power.size());
		
		for(size_t col = 0; col < power.size(); col++)
			//if the indices want to set the lower triangle of the matrix, the powers will be set to 0 to avoid misscalculations
			if(col < row)
			{
				vector<int> tmp(power[0].size(), 0);
				templatemat[row][col] = tmp;
			}
			else
				//if in upper triangle, add the powers as a result of the superposition of the two jacobian matrices
				templatemat[row][col] = add_powers(power[col], power[row]);
			
	}
	
	//if done, check for the highest occuring powers
	set_highest_powers();
}

vector<int> CovMatHandler::add_powers(const vector<int>& a, const vector<int>& b)
{
	vector<int> tmp;
	
	for(size_t i = 0; i < a.size(); i++)
		tmp.push_back(a[i] + b[i]);
	
	return tmp;
	
}

vector<vector<double>> CovMatHandler::get_specific_mat(vector<double>& parameters)
{
	vector<vector<double>> powers = get_powers(parameters);
	
	vector<vector<double>> result;
	
	result.resize(templatemat.size());
	const size_t npars = powers.size();		
	if(!is_diag)
		for(size_t row = 0; row < templatemat.size(); row++)
		{
			result[row].resize(templatemat[row].size());
			
			for(size_t col = 0; col < templatemat[row].size(); col++)
				result[row][col] = get_specific_element(row, col, powers, npars);

		}
	else
		for(size_t i = 0; i < templatemat.size(); i++)
		{
			result[i].resize(templatemat[i].size());
			result[i][i] = get_specific_element(i, i, powers, npars);
		}
	
	return result;
}

void CovMatHandler::get_specific_mat(vector<double>& parameters, vector<vector<double>>& result)
{

	vector<vector<double>> powers = get_powers(parameters);
	
	
	if(result.empty())
	{
		result.resize(templatemat.size());
		
		for(size_t row = 0; row < templatemat.size(); row++)
			result[row].resize(templatemat[row].size());
	}
	
	const size_t templatematsize = templatemat.size();
	const size_t npars = powers.size();	
	if(!is_diag)
		for(size_t row = 0; row < templatematsize; row++)			
			for(size_t col = row; col < templatematsize; col++)
				result[row][col] = get_specific_element(row, col, powers, npars);
					
	else
		for(size_t i = 0; i < templatematsize; i++)
			result[i][i] = get_specific_element(i, i, powers, npars);
		
}

inline double CovMatHandler::get_specific_element(size_t row, size_t col, const vector<vector<double>>& powers, const size_t npars)
{
	double element = 1;
	


	for(size_t i = 0; i < npars; i++)
		element *= powers[i][templatemat[row][col][i]];

	
	return element;
}
	
double CovMatHandler::get_var(const vector<vector<double>>& mat, size_t index)
{
	double result = 0;
	
	const size_t nrows = covmat[index].size();
	size_t row, col;
	
	if(!is_diag)
		for(row = 0; row < nrows; row++)
			for(col = row; col < nrows; col++)
				result += mat[row][col] * covmat[index][row][col];
				
	else
		for(size_t i = 0; i < nrows; i++)
			//~ result += mat[i][i] * mat[i][i] * covmat[index][i][i] * covmat[index][i][i]; //editiert fuer gauss
			result += mat[i][i] * covmat[index][i][i];
	
	return result * result; //modified
}

void CovMatHandler::set_highest_powers()
{
	highest_powers.assign(templatemat[0][0].size(), 0);
	
	for(size_t row = 0; row < templatemat.size(); row++)
		for(size_t col = 0; col < templatemat[row].size(); col++)
			for(size_t i = 0; i < highest_powers.size(); i++)
				if(templatemat[row][col][i] > highest_powers[i])
					highest_powers[i] = templatemat[row][col][i];
	
}

vector<vector<double>> CovMatHandler::get_powers(vector<double>& parameters)
{
	vector<vector<double>> powers;
	
	powers.resize(highest_powers.size());
	for(size_t i = 0; i < powers.size(); i++)
		powers[i].resize(((size_t) highest_powers[i]) + 1);
		

	for(size_t i = 0; i < powers.size(); i++)
		for(size_t j = 0; j < powers[i].size(); j++)
			if(j == 0)
				powers[i][0] = 1;
			else
				powers[i][j] = powers[i][j - 1] * parameters[i];
	
	return powers;
}

void CovMatHandler::set_diagmat(vector<double> vec, size_t index)
{
	covmat[index].resize(vec.size());
	
	for(size_t i = 0; i < vec.size(); i++)
	{
		covmat[index][i].assign(vec.size(), 0);
		covmat[index][i][i] = vec[i];
	}
	
	is_diag = 1;
	
}

void CovMatHandler::is_symm(size_t index){
	
	for(size_t i = 0; i < covmat[index].size(); i++)
		for(size_t j = 0; j < covmat[index][i].size(); j++)
			if(fabs(covmat[index][i][j] / covmat[index][j][i] - 1) > 0.01)
			{
				cout << "asymmetry in " << index << " " << covmat[index][i][j] / covmat[index][j][i] << endl;
				return;
			}
	
	
}	
	
void CovMatHandler::make_triangle(size_t index){
	
	for(size_t row = 0; row < covmat[index].size(); row++)
		for(size_t col = 0; col < row; col++)
			covmat[index][col][row] += covmat[index][row][col];
			
}
