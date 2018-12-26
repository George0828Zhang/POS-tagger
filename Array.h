#include <vector>
#include <numeric>
#include <cassert>

using namespace std;
template <class T>
class
Array
{
	vector<T> _data;
	int _nelem;
	vector<int> _access;
	vector<int> _shape;

	public:
		Array()
			: _nelem(0)
		{
			
		}
		Array(const std::vector<int>& shape)
			: _nelem(std::accumulate(shape.begin(), shape.end(),
					1, multiplies<int>())),
			_shape(shape)
		{
			_data.resize(_nelem, 0);
			_access.resize(shape.size(), 0);
			int val = 1;
			for(int i = shape.size() - 1; i >= 0; i--){				
				_access[i] = val;
				val *= shape.at(i);
			}
		}

		~Array()
		{
		}

		T& operator[](const vector<int>& indices)
		{
			int ptr = 0;
			for(int i = 0; i < _access.size(); i++)
				ptr += _access[i]*indices[i];
			
			return _data[ptr];
			// return _data.at(ptr);
			// return _data.at(inner_product(_access.begin(), _access.end(), indices.begin(), 0));
		}

		void reshape(const vector<int>& new_shape)
		{
			int newsz = accumulate(new_shape.begin(), new_shape.end(), 1, multiplies<int>());
			if(newsz > _data.size()){
				_data.resize(newsz, 0);
			}
			_nelem = newsz;

			int newsz2 = new_shape.size();	
			_shape = vector<int>(new_shape);
			_access.resize(newsz2, 0);
			// if(newsz > _access.size()){
			// 	_access.resize(newsz, 0);
			// }
			int val = 1;
			for(int i = newsz2 - 1; i >= 0; i--){				
				_access[i] = val;
				val *= new_shape.at(i);
			}
		}

		int size() const
		{
			return _nelem;
		}
};
