#include <vector>
#include <numeric>

using namespace std;
template <class T>
class
Array
{
	vector<T> _data;
	int _nelem;
	vector<int> _access;

	public:
		Array()
			: _nelem(0)
		{
			
		}
		Array(const std::vector<int>& shape)
			: _nelem(std::accumulate(shape.begin(), shape.end(),
					1, multiplies<int>()))
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
			return _data.at(inner_product(_access.begin(), _access.end(), indices.begin(), 0));
		}

		void reshape(const vector<int>& new_shape)
		{
			int newsz = accumulate(new_shape.begin(), new_shape.end(), 1, multiplies<int>());
			if(newsz > _data.size()){
				_data.resize(newsz, 0);
			}
			_nelem = newsz;

			newsz = new_shape.size();
			if(newsz > _access.size()){
				_access.resize(newsz, 0);
			}
			int val = 1;
			for(int i = newsz - 1; i >= 0; i--){				
				_access[i] = val;
				val *= new_shape.at(i);
			}
		}

		int size() const
		{
			return _nelem;
		}
};