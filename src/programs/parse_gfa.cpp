
#include <tuple>

#include <sstream>
#include <fstream>
#include <vector>
#include <set>
#include <deque>
#include <tuple>
#include <map>

#include "../lib/pos_EDS_to_GFA.cpp"

namespace gedmap_parse_gfa{
	using namespace std;
	using namespace sdsl;

	vector<string> seq_names;
	uint32_t ORIG_NODE_NUMBER;
	uint64_t MAX_OFFSET;

	template<class T>
	void append_to(vector<T> & A, vector<T>  & B){
		A.insert( A.end(), B.begin(), B.end());
	};

	struct edge_t{
		uint32_t target;
		bool from_orient;
		bool to_orient;

		edge_t(){};
		edge_t(uint32_t target,bool from_orient,bool to_orient):target(target),from_orient(from_orient),to_orient(to_orient){};
		edge_t(uint32_t target,char from_orient_c,char to_orient_c):target(target){

			from_orient = (from_orient_c == '+')?PLUS:!PLUS;
			to_orient   = (to_orient_c   == '+')?PLUS:!PLUS;
		};

		bool equal( edge_t & e2){
			return target == e2.target && from_orient == e2.from_orient && to_orient == e2.to_orient;
		}
		void set_vals(edge_t & e2){
			target 		= e2.target ;
			from_orient	= e2.from_orient;
			to_orient		= e2.to_orient;
		}

		static const bool PLUS = true;
	};

	struct node_t{
		string seq;
		char type;

		string ids;

		//FOR TRACING BACK TO THE NODE
		vector<uint32_t>	orig_seq_idx;
		vector<uint64_t>	orig_seq_pos;
		vector<string>		orig_nodes;
		vector<bool> 		orientation;
		vector<bool>		node_start;

		node_t(){type = 0;};

		//node_t(string seq, string ids):seq(seq), ids(ids){type = ALIVE;};

		node_t(string seq, string id, uint32_t name_idx, uint64_t pos):seq(seq){
			type = ALIVE;
			ids = id;
			orig_seq_idx = vector<uint32_t>(1,name_idx);
			orig_seq_pos = vector<uint64_t>(1,pos);
			orig_nodes   = vector<string>(1,id);
			orientation  = vector<bool>(1,false);
			node_start   = vector<bool>    (seq.size(),0);
			node_start[0] = 1;
			if(pos > MAX_OFFSET) MAX_OFFSET = pos;
		};

		static const char FW = 1;
		static const char BW = 2;
		static const char FUSED = 4;
		static const char ALIVE = 8;

		bool alive(){
			return type & ALIVE;
		}
		void kill(){
			if( alive() ) type -= ALIVE;
			if( alive() ) throw std::runtime_error("KILL: ALREADY DEAD!");
		}

		void reverse(){
			if(orig_nodes.size() != 1)
			throw runtime_error("cannot reverse node that has not one originating node");

			seq = gedmap_encode::rev_complement(seq);
			if      (type & node_t::FW && !(type & node_t::BW))	type = type - node_t::FW + node_t::BW;
			else if (type & node_t::BW && !(type & node_t::FW))	type = type - node_t::BW + node_t::FW;

			ids =  '*' + ids;
			orientation[0] = !orientation[0];
			// -> orig_nodes and node_start do not change
		}

		node_t copy(){
			node_t c;
			c.seq	= seq;
			c.type	= type;
			c.ids	= ids;
			c.orig_seq_idx = orig_seq_idx;
			c.orig_seq_pos = orig_seq_pos;
			c.orig_nodes	= orig_nodes;
			c.orientation	= orientation;
			c.node_start	= node_start;
			return c;
		}

		node_t reverse_copy(){
			node_t c = copy();
			c.reverse();
			return c;
		}

		void append(node_t & v){
			seq += v.seq;
			ids += v.ids;
			append_to(orig_seq_idx, v.orig_seq_idx);
			append_to(orig_seq_pos, v.orig_seq_pos);
			append_to(orig_nodes,   v.orig_nodes);
			append_to(orientation,  v.orientation);
			append_to(node_start,   v.node_start);
		}
	};

	using adj_t = vector< pair< vector<edge_t>, vector<edge_t>>>;

	tuple<uint32_t,uint32_t,uint32_t> count_nodes_edges(string fname_gfa){
		ifstream gfa_in = ifstream(fname_gfa, ifstream::in);
		std::string line;
		uint32_t nc = 0, ec = 0;
		uint32_t lines = 0;
		while(std::getline(gfa_in, line)){
			if(line.size()>0 && line[0]=='S') nc++;
			if(line.size()>0 && line[0]=='L') ec++;
			lines++;
			if(!(lines%23957)) cout << "\r" << lines << " lines" << flush;
		}
		gfa_in.close();
		cout << "\r";

		ORIG_NODE_NUMBER = nc;
		return make_tuple(nc,ec,lines);
	}

	void all_about(size_t i, vector<node_t> & V, adj_t & A, bool init=true){

		cout << "NODE " << i << endl;
		cout << "-seq= " << V[i].seq << endl;
		cout << "-type= " << (int) V[i].type<< endl;
		cout << "-ids= " << V[i].ids << endl;
		cout << "INCOMING" << endl;
		for(auto e : A[i].first)  cout << "<" << e.target << " " << e.from_orient << "->" << e.to_orient << endl;
		cout << "OUTGOING" << endl;
		for(auto e : A[i].second) cout << ">" << e.target << " " << e.from_orient << "->" << e.to_orient << endl;
		if(!init) return;
		cout << "ADJ_NODES" << endl;
		for(auto e : A[i].first)  all_about(e.target,V,A,false);
		for(auto e : A[i].second) all_about(e.target,V,A,false);

	}

	void print_graph(vector<node_t> & V, adj_t & A){
		for(size_t u = 0; u < V.size(); u++ )
			all_about(u,V,A,false);
	}

	void stats(vector<node_t> & V, adj_t & A){
		cout << "##################### GRAPH STATS  ###############" << endl;
		cout << "- node count                        " << V.size() << endl;

		uint32_t fw_nodes = 0, bw_nodes = 0, fwbw_nodes = 0, no_edge_nodes = 0, dead_nodes = 0, fused_nodes = 0 , sources = 0, sinks = 0;

		for( auto u : V){
			if(!u.alive()){ dead_nodes++; continue;}
			bool bw = u.type & node_t::BW;
			bool fw = u.type & node_t::FW;

			if(  bw &&  fw) fwbw_nodes++;
			if(  bw && !fw) bw_nodes++;
			if( !bw &&  fw) fw_nodes++;
			if( !bw && !fw) no_edge_nodes++;

			if( u.type & node_t::FUSED) fused_nodes++;

		}

		cout << "- nodes only in forward   direction " << fw_nodes << endl;
		cout << "- nodes only in backwards direction " << bw_nodes << endl;
		cout << "- nodes used in both directions     " << fwbw_nodes << endl;
		cout << "- nodes without edges               " << no_edge_nodes << endl;
		cout << "- nodes dead                        " << dead_nodes << endl;
		cout << "- nodes with fusions                " << fused_nodes << endl;

		uint32_t edge_count_in  = 0;
		uint32_t edge_count_out = 0;
		for( size_t i = 0 ; i < A.size(); i++ )
			if(V[i].alive()){
				edge_count_in  += A[i].first .size();
				edge_count_out += A[i].second.size();
				if(A[i].first.empty()) sources++;
				if(A[i].second.empty()) sinks++;
			}

		cout << "- sources                           " << sources << endl;
		cout << "- sinks                             " << sinks << endl;
		cout << "- edges count in                    " << edge_count_in << endl;
		cout << "- edges count out                   " << edge_count_out << endl;
		cout << "##################################################" << endl;

		for(size_t i = 0; i < V.size(); i++){
			if( ! V[i].alive() ) continue;

			for( auto e : A[i].second )
				if ( ! V[e.target].alive() )
				{
					cout << "ERROR alive->dead " << i << " -> " << e.target <<  " " << V[i].ids << " -> " << V[e.target].ids << endl;
					all_about(i,V,A);
					all_about(e.target,V,A);
					exit(0);
				}

		}

	}

	void replace_edge( vector<edge_t> & E , edge_t  old_e, edge_t  new_e ){
		for( auto & e : E )
			if(e.equal(old_e))
				e.set_vals(new_e);
	}

	bool bubble_between(size_t u, size_t v, vector<node_t> & V , adj_t & A, uint32_t width, uint32_t max_path_count){
		set   <size_t> I; //inner nodes
		vector< pair<size_t,uint32_t> > S;
		uint32_t observed_paths = 0;

		S.emplace_back( u , 0);

		while( ! S.empty() ){
			size_t e; uint32_t s; tie(e,s) = S.back(); S.pop_back();
			if( e != v){
				auto & suc = A[e].second;
				if( suc.empty()) return false;
				if( suc.size() + observed_paths + S.size() > max_path_count) return false;

				if( e!=u ) {
					if(V[e].type & node_t::FUSED) return false;
					else I.insert(e);
				}
				for(auto x : suc) {
					if(x.target == e) return false;
					if( s + V[x.target].seq.size() > width)   return false;
					S.emplace_back(x.target, s + V[x.target].seq.size());
				}
			}else{
				observed_paths++;
			}
		}

		//check_if_bubbly
		uint32_t sum_of_outgoing_edges = A[u].second.size() ;
		uint32_t sum_of_incoming_edges = A[v].first .size() ;
		for( auto x : I ) sum_of_outgoing_edges += A[x].second.size();
		for( auto x : I ) sum_of_incoming_edges += A[x].first .size();

		if(sum_of_incoming_edges != sum_of_outgoing_edges) return false;

		return true;
	}

	vector< vector <size_t> > paths_between(size_t u, size_t v,  adj_t & A){
		deque < vector<size_t> > Q;
		vector< vector<size_t> > P;

		for( auto e : A[u].second)
			if(e.target == v)	P.push_back( vector<size_t>(0) );
			else				Q.push_back( vector<size_t>(1,e.target) );

		while( !Q.empty() ){
			size_t last = Q.front().back();

			if(last == v){ //BASE CASE, PATH REACHED v
				P.emplace_back(  Q.front() );
				P.back().pop_back();
				Q.pop_front();
				continue;
			}

			if( A[last].second.size() == 1 ) //SPECIAL CASE JUST TO BE MORE EFFICIENT
				Q.front().emplace_back( A[last].second.front().target );
			else{
				for( auto e : A[last].second ){
					Q.emplace_back ( Q.front());
					Q.back().emplace_back(e.target);
				}
				Q.pop_front();
			}
		}
		return P;
	}

	char classify_bubble(vector< vector <size_t> >  & P,  vector<node_t> & V){
		bool all_one     = true;
		bool all_max_one = true;
		bool snp         = true;

		if( P.size() == 1 && P.front().size() == 0) return 'C';
		for( auto p : P){
			all_one     = all_one     && (p.size() == 1);
			all_max_one = all_max_one && (p.size() <= 1);
			snp         = snp && all_one && (V[ p.front() ].seq.size() == 1);
		}

		if(snp)     return 'S';
		if(all_one) return 'A';
		if(all_max_one && P.size() == 2) return 'I';
		return 'B';

	}

	node_t bubble_to_brackets(vector< vector <size_t> >  & P,  vector<node_t> & V){
		size_t node_count = 0;
		for(auto p: P) node_count += p.size();
		if(node_count == 0) return node_t();

		node_t tmp;

		for(auto p : P){ //ITERATE OVER WALKS
			tmp.seq += '|'; tmp.ids += '|'; tmp.node_start.push_back(0);
			for(size_t i = 0; i < p.size(); i++){
				tmp.append( V[p[i]]);
				if( i < p.size()-1 && i > 0 ) tmp.ids += ',';
			}
		}
		tmp.seq[0] = '('; tmp.ids[0] = '(';
		tmp.seq   += ')'; tmp.ids   += ')'; tmp.node_start.push_back(0);
		return tmp;
	}

	char bubble_start( size_t u,  vector<node_t> & V , adj_t & A, uint32_t max_width, uint32_t max_path_count){
		size_t last = u;
		size_t path_size = 0;

		while(true){
			auto & suc = A[last].second;
			if( suc.empty() ) return 0; //PATH ENDS BEFORE BUBBLE END
			last = suc.front().target;
			if(last == u) return 0;
			if( bubble_between(u,last,V,A,max_width,max_path_count) ) break;

			path_size += V[last].seq.size();
			if( path_size > max_width) return 0; //PATH IS TO BIG
		}

		vector< vector <size_t> >  P = paths_between(u, last, A);

		for( auto p : P) for (auto x : p) if(V[x].type & node_t::FUSED) return 0;

		char TYPE  = classify_bubble(P,V);

		//EDIT GRAPH
		if(TYPE != 'C'){
			V[u].type |= node_t::FUSED;
			node_t tmp = bubble_to_brackets(P,V);
			V[u].append(tmp);
			for( auto p : P) for( auto x : p) V[x].kill();
		}else{
			V[u].ids += ',';
			for( auto p : P) for( auto x : p) if(V[x].type & node_t::FUSED) V[u].type |= node_t::FUSED;
			if(V[last].type & node_t::FUSED) V[u].type |= node_t::FUSED;
		}

		V[u].append(V[last]);
		V[last].kill();

		A[u].second = move(A[last].second);
		for( auto e : A[u].second)
			replace_edge( A[e.target].first, edge_t( last, true, true )  , edge_t( u, true, true )  );

		return TYPE;
	}

	void adjust_edges_after_split( size_t u, size_t v,  adj_t & A){

		vector<edge_t> u_in;
		vector<edge_t> u_out;
		vector<edge_t> v_in;
		vector<edge_t> v_out;

		//SPLIT EDGE VECTOR IN TWO VECTORS

		for(auto e : A[u].first) //ITERATE INCOMING EDGES
			if(e.to_orient == edge_t::PLUS)	u_in.push_back(e);
			else							v_in.push_back(e);


		for(auto e : A[u].second) //ITERATE OUTGOING EDGES
			if(e.from_orient == edge_t::PLUS)	u_out.push_back(e);
			else							v_out.push_back(e);

		//ADJUST DIRECTION OF SECOND VECTOR

		for(auto &e : v_in)  e.to_orient   = !e.to_orient;
		for(auto &e : v_out) e.from_orient = !e.from_orient;


		u_in .shrink_to_fit();
		u_out.shrink_to_fit();
		v_in .shrink_to_fit();
		v_out.shrink_to_fit();

		//ADJUST EDGES IN OTHER NODES

		for(auto e : v_in)  replace_edge( A[e.target].second , edge_t( u ,  e.from_orient, !e.to_orient), edge_t( v , e.from_orient, e.to_orient)  );
		for(auto e : v_out) replace_edge( A[e.target].first  , edge_t( u , !e.from_orient,  e.to_orient), edge_t( v , e.from_orient, e.to_orient)  );



		A[u].first .swap(u_in);
		A[u].second.swap(u_out);
		A[v].first .swap(v_in);
		A[v].second.swap(v_out);

	}

	void add_nodes_for_rev_compl(vector<node_t> & V , adj_t & A){
		size_t count_old = V.size();
		size_t count_new = 0;
		for(auto u : V) if( u.type & node_t::FW && u.type & node_t::BW) count_new++;
		//make space;
		V.resize(count_old+count_new);
		A.resize(count_old+count_new);
		//
		cout << "add " << count_new << " nodes for reverse compliments" << endl;
		//
		size_t vp = count_old;
		for(size_t up = 0 ; up < count_old; up++) {
			auto u = &(V[up]);
			if( u->type & node_t::FW && u->type & node_t::BW){
				u->type -= node_t::BW;
				V[vp] = u->reverse_copy();
				adjust_edges_after_split( up , vp , A);
				vp++;
			}
		}
	}

	/** @brief reorientate edges of a node */
	void reorientate_edges_of_node(size_t i, adj_t & ADJ){
		//INCOMING EDES
		for( auto& e : ADJ[i].first ){
			e.to_orient = !e.to_orient;
			//FIND COPY OF THIS EDGE
			for(auto& e2 : ADJ[e.target].second )
				if(e2.target == i)
					e2.to_orient = !e2.to_orient;
		}
		//OUTGOING EDES
		for( auto& e : ADJ[i].second ){
			e.from_orient = !e.from_orient;
			//FIND COPY OF THIS EDGE
			for(auto& e2 : ADJ[e.target].first )
				if(e2.target == i)
					e2.from_orient = !e2.from_orient;
		}
 	}

 	/** @brief nodes that are only used in backward direction are reorientet (label gets its reverse compliment)  */
	void reorientate_nodes(vector<node_t> & V , adj_t & ADJ){
		uint32_t count = 0;
		for(size_t i = 0 ; i < V.size(); i++){
			if( (V[i].type & node_t::BW) &&  !(V[i].type & node_t::FW) ){
				count++;
				reorientate_edges_of_node(i,ADJ);
				V[i].reverse();
			}
		}
		cout << count << " nodes reoriented" << endl;
	}

	pair<string,uint64_t> get_seq_data_from_line( std::istringstream & iss){
		string tmp, name = ""; uint64_t off = 0;
		while(iss){
			iss >> tmp;
			if(tmp.rfind("SN:Z:", 0) == 0) name =       tmp.substr(5);
			if(tmp.rfind("SO:i:", 0) == 0) off  = stoul(tmp.substr(5));
		}
		return make_pair(name,off);
	};

	template<typename T>
	size_t get_index( vector<T> & v, T t, size_t hint = 0 ){
		if(hint < v.size() && v[hint] == t) return hint;
		for(size_t i = 0; i < v.size(); i++) if(v[i] == t) return i;
		v.push_back(t); return v.size()-1;
	}

	/** @brief reads GFA file */
	pair<vector<node_t>, adj_t > read_gfa( string fname_gfa){

		seq_names = vector<string>(0);
		size_t seq_names_idx = -1;

		std::string line;

		cout << "scan file: " << flush;
		uint32_t node_count,edge_count,line_count, nc = 0, ec = 0, lines = 0;
		tie(node_count, edge_count,line_count) = count_nodes_edges(fname_gfa);
		cout << node_count << "nodes and "<< edge_count << " edges" << endl;
		size_t tot_size = 0;

		vector<node_t> V;
		V.reserve(node_count+1);
		map<string,size_t> name_to_Vid;

		adj_t ADJ ( node_count+1  );

		ifstream GFA = ifstream(fname_gfa, ifstream::in);

		cout << "read nodes and edges       " << endl;
		cout << std::setprecision(2) << std::fixed;

		while(std::getline(GFA, line)) {     // '\n' is the default delimiter
			lines++;
			if(!(lines%23957)) cout << "\r" << (nc/(double)(node_count))*100 << "% nodes " << (ec/(double)(edge_count))*100 <<  "% edges " <<flush;

			std::istringstream iss(line);

			char T;
			iss >> T;

			switch (T){
				case '#': //COMMENT
					break;
				case 'H': //HEADER
					break;
				case 'S': {//SEGMENT
					nc++;
					std::string name,seq;
					iss >> name >> seq;
					tot_size += seq.size();

					string seq_name; uint64_t seq_off;
					tie(seq_name, seq_off) = get_seq_data_from_line(iss);
					seq_names_idx = get_index(seq_names, seq_name,  seq_names_idx);

					if(name_to_Vid.emplace(name,V.size()).second) //insert key-val pair to map (if false, key already exists)
						V.emplace_back(move(seq),name , seq_names_idx, seq_off);
					else
						throw invalid_argument("Multiple segments with name >" + name + "< in gfa file");
					break;
				}
				case 'L': {//LINK
					ec++;
					std::string from_s, to_s, overlap;
					char from_o, to_o;

					iss >> from_s >> from_o >> to_s >> to_o >> overlap;
					size_t from, to;
					try{
						from = name_to_Vid.at(from_s);
						to   = name_to_Vid.at(to_s);
					}catch(std::out_of_range & oor){
						throw invalid_argument("Unspecified segment: Link: " + from_s + " > " + to_s  + " uses unspecified name in From or To field.");
					}

					ADJ[to]  .first.emplace_back(from, from_o, to_o);
					ADJ[from].second.emplace_back(to , from_o, to_o);

					V[from].type |= (from_o=='+'?node_t::FW : node_t::BW);
					V[to]  .type |= (to_o  =='+'?node_t::FW : node_t::BW);

					edge_count++;
					break;
				}
				case 'J': //JUMP
				case 'C': //CONTAINMENT
				case 'P': //PATH
				case 'W': //WALK
				default: break;
			}
		}
		cout << "\r";
		return make_pair(V,ADJ);
	}

	void transform_to_eds_graph(vector<node_t> & V , adj_t & A){

		uint32_t alts = 0, bubbles = 0, snps = 0, indels = 0, fuse = 0 ;

		for(size_t u = 0; u < V.size(); u++){
			if(!(u%23957)) cout << "\r" << (u / (double) V.size())*100 << "%" << flush;

			if( ! V[u].alive() ) continue;
			char t = bubble_start(u,V,A,BUB_max_length, BUB_max_path_c);
			switch (t){
				case 'S': snps++;    break;
				case 'I': indels++;  break;
				case 'C': fuse++;    break;
				case 'A': alts++;    break;
				case 'B': bubbles++; break;
			}
			if(t) u--;
		}
		cout << "\r-snps           " << snps
			<< "\n-bubbles        " << bubbles
			<< "\n-alternatives   " << alts
			<< "\n-indels         " << indels
			<< "\n-concatenations " << fuse << endl;

	}

	void clean_graph(vector<node_t> & V, adj_t & A){
		size_t new_node_count = 0;
		for(auto v : V) if( v.alive() ) new_node_count++;

		vector<size_t> new_index(V.size(), -1);
		vector<node_t> new_V(1); // (1) to leave an empty spot, to avoid a zero index
		adj_t new_A(1);
		new_V.reserve(new_node_count+1);
		new_A.reserve(new_node_count+1);

		if(print) cout << "move nodes" << endl;

		for(size_t u = 0; u < V.size(); u++)
			if( V[u].alive() ){
				new_index[u] = new_V.size();
				new_V.push_back( move(V[u]) );
				new_A.push_back( move(A[u]) );

			}
		if(print) cout << "correct edges " << endl;

		for(size_t u = 0; u < new_V.size(); u++)
			for(auto & e : new_A[u].second){
				e.target = new_index[e.target];
				new_V[e.target].type |= (e.to_orient==edge_t::PLUS)? node_t::FW : node_t::BW;
				new_V[u].type       |= (e.from_orient==edge_t::PLUS)? node_t::FW : node_t::BW;
			}
		V.swap(new_V);
		A.swap(new_A);
	}

	void write_gfa_graph( vector<node_t> & V, adj_t & A, ofstream & off){
		off << "# This is a graph generated by gedmap\n";
		off << "# The forth column of a segment entry displays the number of the corresponding nodes in the ogiginal gfa file.\n";
		off << "H	VN:Z:1.0\n";
		for( size_t i = 1 ; i < V.size(); i++ )
			off << "S\t" << i << '\t' << V[i].seq << '\t' << V[i].ids << '\n';

		for( size_t i = 1 ; i < A.size(); i++ )
			for(auto e : A[i].second)
				off << "L\t" << i << (e.from_orient?"\t+\t":"\t-\t") << e.target << (e.to_orient?"\t+\n":"\t-\n");

		off.flush();
	}

	uint32_t bits8(uint64_t x){
		uint32_t y = bitsneeded(x);
		uint32_t z = 8;
		while(z < y) z+=8;
		return z;
	}

	template<unsigned char X,class T>
	void write_to( sdsl::int_vector<X> & iv , size_t pos, vector<T> & v){
		if(v.size() + pos > iv.size()) throw runtime_error("write_to: cannot copy to iv: ");
		for(auto t : v) iv[ pos++ ] = t;
	}

	void write_EDSG( vector<node_t> & V, adj_t & A, string out){

		ofstream adj_out  (out + "." + FEX_ADJ,  std::ios::binary | std::ios::trunc | std::ios::out);
		ofstream edsg_out (out + "." + FEX_EDS,  std::ios::out);
		ofstream p2gfa_out(out + "." + FEX_2GFA, std::ios::out);

		if(print) cout << "prepare output" << endl;
		size_t seq_lenght = V.size();
		size_t seq_count = 0;
		for( auto v : V) if(v.alive()) seq_lenght += v.seq.size();
		for( auto v : V) seq_count += v.orientation.size();


		sdsl::int_vector<1> node_start_ind(seq_lenght,0,1);
		//ENTRY PER NODE
		vector<string> node_name; node_name.reserve(ORIG_NODE_NUMBER);
		sdsl::int_vector<1> orientation (seq_count, 0, 1);
		sdsl::int_vector<0> seq_name_idx(seq_count, 0, bits8(seq_names.size())); //SEQ IDX IN SEQ_NAME
		sdsl::int_vector<0> offset      (seq_count, 0, bits8(MAX_OFFSET));       //SEQ OFFSET OF NODE



		vector<uint64_t> starting_hash(V.size()+1);

		size_t cur_node = 0;

		if(print) cout << "write GEDS" << endl;

		for( size_t i = 0; i < V.size(); i++){
			if(!V[i].alive()) continue;
			starting_hash[i] = edsg_out.tellp();
			edsg_out << '#';

			write_to(node_start_ind, edsg_out.tellp(), V[i].node_start);

			edsg_out << V[i].seq;

			append_to(node_name, V[i].orig_nodes);
			write_to(orientation,  cur_node, V[i].orientation);
			write_to(seq_name_idx, cur_node, V[i].orig_seq_idx);
			write_to(offset,       cur_node, V[i].orig_seq_pos);

			if(V[i].orig_nodes.size() != V[i].orientation.size()
				||V[i].orig_nodes.size() != V[i].orientation.size()
				||V[i].orig_nodes.size() != V[i].orig_seq_idx.size()
				||V[i].orig_nodes.size() != V[i].orig_seq_pos.size())
				throw runtime_error("sizes do not match");

			cur_node += V[i].orig_nodes.size();
		}

		starting_hash.back() = edsg_out.tellp();
		edsg_out << '#';
		edsg_out.close();


		if(print) cout << "write 2gfa" << endl;

		{
			pos_EDS_to_GFA_type p2gfa(
				move(node_name),
				move(orientation),
				move(seq_name_idx),
				move(offset),
				move(seq_names),
				move(node_start_ind));
			sdsl::serialize( p2gfa,p2gfa_out);
			p2gfa_out.close();
		}

		if(print) cout << "write adj" << endl;

		vector<adjacency::edge_type> Edges;
		for( size_t i = 0; i < A.size(); i++ )
			for( auto e : A[i].second)
				Edges.emplace_back ( starting_hash[i+1] , starting_hash[e.target] );

		adjacency adj (Edges);
		sdsl::serialize( adj,adj_out);
		adj_out.close();




	}

	int main(int argc, char** argv){
		using namespace gedmap_io;

		sys_timer t;
		t.start();

		print_prog_headline("GEDMAP PARSE");
		//OUTPUT INFO
		print_row("Parse GFA: " , argv[2]);
		dotline();

		string fname_gfa;
		string out_prefix;

		handle_input(argc,argv,fname_gfa,out_prefix);


		print_row("read gfa");
		vector<node_t> V ;
		adj_t A;
		tie(V,A) =  read_gfa( fname_gfa);
		stats(V,A);

		print_row("reorientate nodes, that are only used in backwards direction");
		reorientate_nodes(V,A);
		if(print) stats(V,A);

		print_row("split nodes, that are used in both directions");
		add_nodes_for_rev_compl(V,A);
		if(print) stats(V,A);

		print_row("detect bubbles");
		transform_to_eds_graph(V,A);
		if(print) stats(V,A);
		print_row("cleanup");
		clean_graph(V,A);
		if(print) stats(V,A);


		print_row("write gfa graph with eds labels (temporary)");
		ofstream tmp(out_prefix+".transformed.gfa");
		write_gfa_graph(V,A,tmp);



		print_row("write EDS graph for gedmap");
		stats(V,A);
		write_EDSG(V,A,out_prefix);




		//WRITE INFO
		dotline();
		print_row("TIME in total: ", t.stop_and_get() , " s" );
		dotline();

		return 0;
	}
}
