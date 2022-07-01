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

	int n_missing=0;
	

	int ploidy=0;

	T& operator[](unsigned i){
		return data[i];
	}

	bool is_empty() const {
		return data == NULL;
	}

	// int do_gt_sfs(){
		// fprintf(stderr,"\n\nHERE!!!\n\n");
	// };

	void destroy(){
		free(data);
		data=NULL;
	}

	~get_data(){
		free(data);
	}

};
//
// template <typename T> struct data_VCF{
//
	 // fprintf(stderr,"\n->->%s\n",hdr->samples[i]);
	// T *data = NULL;
//
	// int size_e=0;
	// int n=0;
	//
//
//
//
	//
//
	// int ploidy=0;
//
	// T& operator[](unsigned i){
		// return data[i];
	// }
//
	// bool is_empty() const {
		// return data == NULL;
	// }
//
	// int do_gt_sfs(){
		// fprintf(stderr,"\n\nHERE!!!\n\n");
	// };
//
	// // printf("REF:%s ALT:%s\n", rec->d.allele[0], rec->d.allele[1]);
	//
//
//
	// void destroy(){
		// free(data);
		// data=NULL;
	// }
//
	// ~get_data(){
		// free(data);
	// }
//
// };
