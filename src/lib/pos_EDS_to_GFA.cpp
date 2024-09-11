#pragma once
#include "pos_EDS_to_FA.cpp"
using namespace std;
using namespace sdsl;




/**
 * @brief data structure that allowes to map a position of the EDS to the original GFA node and position
 * 
 * 
 * node_start_ind = indicator: node_start_ind[i]=1, iff EDS[i] is the start of a node in the original GFA file
 * 

 */
struct pos_EDS_to_GFA_type: public Transform{
	typedef sdsl::int_vector<>::size_type	size_type;
	

	//ENTRY PER NODE
	vector<string> node_name; //NAME IN GFA FILE
	sdsl::int_vector<1> orientation;
	sdsl::int_vector<0> seq_name_idx; //SEQ IDX IN SEQ_NAME
	sdsl::int_vector<0> offset;       //SEQ OFFSET OF NODE

	//ONCE
	vector<string>    seq_name;

	//ENTRY PER CHAR IN EDSG
	sdsl::bit_vector 		node_start_ind;
	sdsl::rank_support_v5<> 	rs_node_start_ind;
	//sdsl::... ss_node_start_ind; //TODO
	
	
	
	/** @brief default constructor */
	pos_EDS_to_GFA_type(){};

	/** @brief parametrized constructor */
	pos_EDS_to_GFA_type(
		vector<string> node_name,
		sdsl::int_vector<1> orientation,
		sdsl::int_vector<0> seq_name_idx,
		sdsl::int_vector<0> offset,
		vector<string>  seq_name,
		sdsl::bit_vector node_start_ind
	):node_name(node_name),orientation(orientation),seq_name_idx(seq_name_idx),offset(offset),seq_name(seq_name),node_start_ind(node_start_ind){
		sdsl::util::init_support(rs_node_start_ind,&node_start_ind);
	};
	
	/** @brief 
	 * <u,p,o> = (pos),
	 * with pos is in node u at position p (with orientation o)
	 */
	std::tuple<std::string,uint64_t,std::string> operator()(uint64_t EDS_pos) const;
	
	bool empty() const;

	size_type 	serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const;
	void 		load(std::istream& in);
};



/** 
 * @brief return position in FFA
 * @return <c,p,n> ,
 * c sequence / node;
 * p = position
 * n = offset???
 */
tuple<std::string,uint64_t,std::string> pos_EDS_to_GFA_type::operator()(uint64_t EDS_pos) const{
	EDS_pos--;
	
	if(EDS_pos >= node_start_ind.size()) throw runtime_error("TRANSFORM: EDS_POS TO BIG");

	uint64_t p = 0;
	while( !node_start_ind[EDS_pos - p] ){
		p++;
		if(EDS_pos < p) throw runtime_error("TRANSFORM: WEIRD");
	};

	size_t node = rs_node_start_ind(EDS_pos-p);
	if(seq_name_idx.size() <= node) throw runtime_error("TRANSFORM: WEIRD2");
	if(seq_name.size() <= seq_name_idx[node]) throw runtime_error("TRANSFORM: WEIRD3");
	string seq = seq_name[ seq_name_idx[node] ];

	if(seq != "" && !orientation[node])	return make_tuple( seq, offset[node]+p, "" );
	if(seq_name_idx.size() <= node) throw runtime_error("TRANSFORM: WEIRD4");
	if(seq == "" && !orientation[node])                  return make_tuple("Node:"+(node_name[node]), p, "");
	return make_tuple("Node_reverse:"+(node_name[node]), p, "");
}

/**
 * RETURNS A FLAT ARRAY THAT CONTAINS THE VALUES OF v AND 0 AS BOARDERS (ALSO IN FIRST AN LAST CELL)
 */
sdsl::int_vector<8> string_vec_to_int_vec(const std::vector<std::string> & v){
	size_t total_len = v.size();
	for( auto &s : v) total_len += s.size();
	sdsl::int_vector<8> out(total_len+1,0,8);
	size_t i = 1;
	for( auto &s : v){
		for(char c : s){
			if( c == 0) throw invalid_argument("string_vec_to_int_vec: string contains \0");
			out[i++] = c;
		}
		i++;
	}
	return out;
}

std::vector<std::string> int_vec_to_string_vec(const sdsl::int_vector<8> & v){

	std::vector<size_t> zero_pos;
	for(size_t i = 0; i < v.size(); i++) if(v[i]==0) zero_pos.push_back(i);
	std::vector<std::string> out(zero_pos.size()-1);
	for(size_t i = 1; i < zero_pos.size(); i++){
		size_t s_start = zero_pos[i-1]+1;
		size_t s_end   = zero_pos[i];
		string s(s_end-s_start,0);
		for(size_t j = 0; j < s.size(); j++)
			s[j] = v[s_start+j];
		out[i-1] = move(s);
	}
	return out;
}




/**
 * @brief serialize pos_EDS_to_FA_type with sdsl methods
 */
pos_EDS_to_GFA_type::size_type pos_EDS_to_GFA_type::serialize(ostream& out, structure_tree_node* v,string name) const{
	
	
	structure_tree_node *child 	= structure_tree::add_child(v, name, sdsl::util::class_name(*this));
	
	pos_EDS_to_GFA_type::size_type written_bytes 	= 0;

	written_bytes += string_vec_to_int_vec(node_name).serialize(out, child, "node_name");
	written_bytes += orientation	.serialize(out, child, "orientation");
	written_bytes += seq_name_idx	.serialize(out, child, "seq_name_idx");
	written_bytes += offset	.serialize(out, child, "offset");
	written_bytes += node_start_ind   .serialize(out, child, "node_start_ind");
	written_bytes += rs_node_start_ind.serialize(out, child, "rs_node_start_ind");

	written_bytes += string_vec_to_int_vec(seq_name).serialize(out, child, "seq_names");

	structure_tree::add_size(child, written_bytes);
	return written_bytes;
}


bool pos_EDS_to_GFA_type::empty() const{
	return (node_name.size() == 0);
}
/**
 * @brief loads pos_EDS_to_FA_type with sdsl methods
 */
void pos_EDS_to_GFA_type::load(std::istream& in){
	sdsl::int_vector<8> tmp;
	tmp				.load(in);
	node_name 		= int_vec_to_string_vec(tmp);
	orientation		.load(in);
	seq_name_idx		.load(in);
	offset			.load(in);
	node_start_ind		.load(in);
	rs_node_start_ind	.load(in,&node_start_ind);
	tmp.load(in);
	seq_name 			= int_vec_to_string_vec(tmp);
}



std::ostream& operator<< (std::ostream& os, const pos_EDS_to_GFA_type& T){
	os << "pos_EDS_to_GFA_type:" << endl;
	os << "node_name ";
	for(auto s : T.node_name) os <<  " " << s;
	os <<  endl;
	os << "orientation " << T.orientation << endl;
	os << "seq_name_idx " << T.seq_name_idx << endl;
	os << "offset " << T.offset << endl;
	os << "node_start_ind " << T.node_start_ind << endl;
	os << "seqs: " << flush;
	for(auto s : T.seq_name) os << "\""<< s << "\""<< endl;
	os <<  endl;
	return os;
}
