/*
 * @template struct
 * @abstract	wrapper for bcf_get_data_*
 *
 * @field data
 * @field size_e	watermark for number of elements
 * @field n			number of returned values
 * 					if <0; error
 * 					else; number of written values
 * 					used for accessing entries
 *
 */
template<typename T> struct get_data{

	T *data = NULL;

	int size_e=0;
	int n=0;

	int n_missing_ind=0;
	

	int ploidy=0;

	T& operator[](unsigned i){
		return data[i];
	}

	bool is_empty() const {
		return data == NULL;
	}

	void destroy(){
		free(data);
		data=NULL;
	}

	~get_data(){
		free(data);
	}

};

struct Strata{
	int nPairs=0;
	char *pairs=NULL;

	int nInds=0;
	char *inds[10];
	int buf_inds=10;
	
	char *id=NULL;
};

struct METADATA{
	Strata S[4];
	int nStrata=0;
	int buf_strata=4;
};
