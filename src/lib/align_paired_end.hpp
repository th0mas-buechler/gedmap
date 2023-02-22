#pragma once

#include "paired_end/edsm_frontier.hpp"

namespace gedmap_align_min {

using pe_event_type = std::tuple<int64_t, uint32_t, std::optional<int64_t>>;
constexpr int64_t pairing_ali_weigth = 10; // ali. dist is multiplied with pairing_ali_weigth for mapping

template< typename int_type, typename MSQ, typename EDGview, typename F, typename GET_TMP_ALIGNMENTS, typename MK_EVENT_POS>
inline
std::vector< temp_alignment<int_type> >
try_mate(
	const MSQ& msq,
	const EDGview& edg_view,
	std::vector< pe_event_type >& events,
	std::vector< pe_event_type >& events_rev,
	std::vector<uint32_t> (&left_mate)[2],
	const std::vector< temp_alignment<int_type> >& tmp_alignments_l, size_t mam,
	F mk_hotspot,
	GET_TMP_ALIGNMENTS get_tmp_alignments,
	MK_EVENT_POS mk_event_pos
) {
#ifndef NDEBUG
	const auto events_old = events;
	const auto events_rev_old = events_rev;
#endif
	const size_t events_pre_size = events.size(), events_rev_pre_size = events_rev.size();
	for (size_t ali_idx = 0; ali_idx < std::min(mam, tmp_alignments_l.size()); ali_idx++) {
		const auto& ali = tmp_alignments_l[ali_idx];
		//tmp_alignments_l[i].source = std::numeric_limits<uint32_t>::max();
		// TODO: find the rightmost matched symbol
		(ali.reverse_compl ? events_rev : events)
			.emplace_back(
				mk_event_pos(ali.eds_pos - ali.read_pos),
				ali_idx,
				pairing_ali_weigth * (int64_t)ali.get_dist());
	}

	// rank candidates according to distance to already aligned "tmp_alignments_l"
	std::vector<std::tuple<uint32_t,bool,double>> opt_mate // (mate index in hotspots_* / pos_pairs_*, r_c, dist)
		( tmp_alignments_l.size()
		, std::make_tuple(std::numeric_limits<uint32_t>::max(), false, std::numeric_limits<double>::infinity()));
	paired_end<int64_t, int64_t>([&opt_mate](auto a, auto b, auto val) {
			if (val < std::get<2>(opt_mate[b]))
				opt_mate[b] = std::make_tuple(a, false, val);
		},
		events,
		edg_view,
		msq,
		PE_FRAGMENT_LENGTH);
	paired_end<int64_t, int64_t>([&opt_mate](auto a, auto b, auto val) {
			if (val < std::get<2>(opt_mate[b]))
				opt_mate[b] = std::make_tuple(a, true, val);
		},
		events_rev,
		edg_view,
		msq,
		PE_FRAGMENT_LENGTH);

	events.resize(events_pre_size), events_rev.resize(events_rev_pre_size);
	assert(events == events_old);
	assert(events_rev == events_rev_old);
	
	for (size_t i = 0; i < 2; i++) {
		left_mate[i].clear();
		left_mate[i].reserve(opt_mate.size());
	}
	return [&] {
		std::vector<hotspot<int_type>> final_hotspots, final_hotspots_r_c;
		for (size_t i_l = 0; i_l < opt_mate.size(); i_l++) {
			const auto&[i_r, r_c, val] = opt_mate[i_l];
			if (val == std::numeric_limits<double>::infinity()) 
				continue;
			left_mate[r_c].emplace_back(i_l);

			(r_c ? final_hotspots_r_c : final_hotspots)
				.emplace_back(mk_hotspot(r_c, i_r, val));
		}
		std::sort(final_hotspots.begin(), final_hotspots.end(),
				[&](const auto& lhs, const auto& rhs) {
					return lhs.quality > rhs.quality;
				});
		std::sort(final_hotspots_r_c.begin(), final_hotspots_r_c.end(),
				[&](const auto& lhs, const auto& rhs) {
					return lhs.quality > rhs.quality;
				});
		return get_tmp_alignments(std::move(final_hotspots), std::move(final_hotspots_r_c));
	}();
}

// events(_rev) contains events for righthand alignment
// this adds the hotspots for the lefthand alignment to the events
// and pairs them
template<typename int_type, typename Environ>
std::vector<std::array< temp_alignment<int_type>, 2> >
fallback(
	const Environ& env,
	std::vector< pe_event_type >& events,
	std::vector< pe_event_type >& events_rev,
	const std::vector< hotspot<int_type> >& hotspots_l,
	const std::vector< hotspot<int_type> >& hotspots_r,
	const std::vector< query_position<int_type> >& pos_pairs_l_r_c,
	const std::vector< query_position<int_type> >& pos_pairs_r_r_c,
	const fasta_read<int_type>& read_l,
	const fasta_read<int_type>& read_r,
	size_t& num_tried,
	double& pairing_time
) {
	constexpr double INF = std::numeric_limits<double>::infinity();
	
	const auto& msq = eget<ED_Graph<uint32_t>>(env);

	const auto& read_l_r_c = read_l.get_rev_compl();
	const auto& read_r_r_c = read_r.get_rev_compl();

	std::vector< std::tuple<double, uint32_t, uint32_t, bool> > mate_candidates; // (val, idx in hotspots, idx in pos_pairs, r_c)
	{
		const size_t events_pre_size = events.size(), events_rev_pre_size = events_rev.size();

		if (not hotspots_l.empty())
		{
			const auto max_qual = hotspots_l[0].quality;
			for (size_t i = 0; i < hotspots_l.size(); i++)
				if (true or hotspots_l[i].quality * 4 >= max_qual * 2)
					events.emplace_back(
						hotspots_l[i].eds_pos - hotspots_l[i].read_pos,
						i,
						-pairing_ali_weigth * (int64_t)hotspots_l[i].quality);
				else break;
		}
		if (not hotspots_r.empty())
		{
			const auto max_qual = hotspots_r[0].quality;
			for (size_t i = 0; i < hotspots_r.size(); i++)
				if (true or hotspots_r[i].quality * 4 >= max_qual * 2)
					events_rev.emplace_back(
						hotspots_r[i].eds_pos - hotspots_r[i].read_pos,
						i,
						-pairing_ali_weigth * (int64_t)hotspots_r[i].quality);
				else break;
		}

		std::vector<std::tuple<uint32_t,double>> opt_mate_l // (mate index in pos_pairs_r_r_c, r_c, dist)
			( hotspots_l.size()
			, std::make_tuple(std::numeric_limits<uint32_t>::max(), INF));
		std::vector<std::tuple<uint32_t,double>> opt_mate_r // (mate index in pos_pairs_l_r_c, dist)
			( hotspots_r.size()
			, std::make_tuple(std::numeric_limits<uint32_t>::max(), INF));

		sys_timer timer;
		timer.start();
		paired_end<int64_t, int64_t>([&opt_mate_l](auto a, auto b, auto val) {
				if (val < std::get<1>(opt_mate_l[b]))
					opt_mate_l[b] = std::make_tuple(a, val);
			},
			events,
			EDG_view<uint32_t>{&msq},
			msq,
			PE_FRAGMENT_LENGTH);
		paired_end<int64_t, int64_t>([&opt_mate_r](auto a, auto b, auto val) {
				if (val < std::get<1>(opt_mate_r[b]))
					opt_mate_r[b] = std::make_tuple(a, val);
			},
			events_rev,
			EDG_view<uint32_t>{&msq},
			msq,
			PE_FRAGMENT_LENGTH);
		pairing_time += timer.stop_and_get();

		events.resize(events_pre_size), events_rev.resize(events_rev_pre_size);

		mate_candidates.reserve( hotspots_l.size() + hotspots_r.size() );
		for (size_t r_c = false; r_c <= true; r_c++)
		{
			auto& opt_mate = r_c ? opt_mate_r : opt_mate_l;
			for (size_t i = 0; i < opt_mate.size(); i++)
			{
				if (std::get<1>(opt_mate[i]) == INF) continue;
				mate_candidates.emplace_back(
					std::get<1>(opt_mate[i]),
					i,
					std::get<0>(opt_mate[i]),
					r_c);
			}
		}
		std::sort(mate_candidates.begin(), mate_candidates.end(), [](const auto& lhs, const auto& rhs) {
			return std::get<0>(lhs) < std::get<0>(rhs);
		});
	}

	// TODO: try aligning all with dist. 0 befor trying with D

	uint32_t D = MAX_DIST.back();

	std::vector< std::array<temp_alignment<int_type>, 2> > tmp_alignments;

	for (size_t i = 0; i < mate_candidates.size() and i < MAX_ALIGNS_T_FALLBACK and tmp_alignments.size() < MAX_ALIGNS_C[0]; i++) {
		num_tried++;

		const auto&[_, hotspot_idx, pos_pair_idx, r_c] = mate_candidates[i];
		const auto&[read_left, read_right] = r_c
			? std::tie(read_r, read_l_r_c)
			: std::tie(read_l, read_r_r_c);
		const auto& hotspot = (r_c ? hotspots_r : hotspots_l)[hotspot_idx];
		const auto& pos_pair = (r_c ? pos_pairs_l_r_c : pos_pairs_r_r_c)[pos_pair_idx];

		// align left
		const auto dist_left = align<false, read_processor::dist_type, uint32_t>(
			eget<std::string>(env),
			eget<adjacency>(env),
			hotspot.eds_pos,
			read_left.sequence, read_left.qual,
			hotspot.read_pos, D);
		temp_alignment<int_type> tmp_ali_left(dist_left, hotspot.eds_pos, hotspot.read_pos, r_c, i);
		if (tmp_ali_left.get_dist() > D) continue;

		const auto dist_right = align<false, read_processor::dist_type, uint32_t>(
			eget<std::string>(env),
			eget<adjacency>(env),
			pos_pair.eds_pos,
			read_right.sequence, read_right.qual,
			pos_pair.read_pos, D);
		temp_alignment<int_type> tmp_ali_right(dist_right, pos_pair.eds_pos, pos_pair.read_pos, r_c, i);
		if (tmp_ali_right.get_dist() > D) continue;

		D = std::max( tmp_ali_left.get_dist(), tmp_ali_right.get_dist() );

		tmp_alignments.push_back({{std::move(tmp_ali_left), std::move(tmp_ali_right)}});

		if (D < gedmap_align_min::DOUBT_DIST)
			break;
	}
	if (tmp_alignments.empty())
		num_tried = 1;

	return tmp_alignments;
}

template<typename int_type, bool reverse = false, typename Environ>
std::array< std::vector<alignment<int_type>>, 2 >
try_pair(
	const Environ& env,
	const fasta_read<int_type>& read_a,
	const fasta_read<int_type>& read_b,
	const fasta_read<int_type>& read_a_r_c,
	const fasta_read<int_type>& read_b_r_c,
	const std::vector<query_position<int_type>>& pos_pairs_a,
	const std::vector<query_position<int_type>>& pos_pairs_b,
	const std::vector<query_position<int_type>>& pos_pairs_a_r_c,
	const std::vector<query_position<int_type>>& pos_pairs_b_r_c,
	size_t& fallback_tries,
	double& pairing_time
) {
	const auto mk_event_pos = [&] (int64_t pos) -> int64_t {
		if constexpr (reverse)
			return eget<std::string>(env).size() - 1 - pos;
		else
			return pos;
	};

	std::vector< pe_event_type > events, events_rev;

	// build events for mate-pairing
	for (uint32_t i = 0; i < pos_pairs_b_r_c.size(); i++)
		events.emplace_back(
			mk_event_pos(pos_pairs_b_r_c[i].eds_pos - pos_pairs_b_r_c[i].read_pos),
			i,
			std::nullopt);
	for (uint32_t i = 0; i < pos_pairs_a_r_c.size(); i++)
		events_rev.emplace_back(
			mk_event_pos(pos_pairs_a_r_c[i].eds_pos - pos_pairs_a_r_c[i].read_pos),
			i,
			std::nullopt);
	
	vector<hotspot<int_type>> hotspots_a, hotspots_b;
	std::vector< temp_alignment<int_type> > tmp_alignments_l, tmp_alignments_r;
	
	size_t last_SPOT_HITS;

		
	std::vector<uint32_t> left_mate[2]; // not-rev, rev

	for (size_t num_try = 0; num_try < std::max(MAX_ALIGNS_T.size(), MAX_ALIGNS_M.size()); num_try++) {
		const size_t align_idx = std::min(num_try, MAX_ALIGNS_T.size() - 1);
		const size_t pair_idx = std::min(num_try, MAX_ALIGNS_M.size() - 1);

		if (num_try == 0 or last_SPOT_HITS != SPOT_HITS[align_idx]) {
			last_SPOT_HITS = SPOT_HITS[align_idx];
			hotspots_a = read_processor::find_hotspots(pos_pairs_a, SPOT_SIZE, last_SPOT_HITS, CHECK_COLLI);
			hotspots_b = read_processor::find_hotspots(pos_pairs_b, SPOT_SIZE, last_SPOT_HITS, CHECK_COLLI);
		}

		read_processor::align_state align_state;
		align_state.D = align_state.D_r_c = MAX_DIST[align_idx];

		tmp_alignments_l = read_processor::align<int_type>//, true>
			( hotspots_a, read_a
			, hotspots_b, read_b
			, env, align_state
			, MAX_ALIGNS_C[align_idx], MAX_ALIGNS_T[align_idx]
			);
		if (tmp_alignments_l.empty())
			continue;

		tmp_alignments_r = try_mate<int_type>(
			[&]() -> std::conditional_t<reverse, EDG_rev<uint32_t>, const ED_Graph<uint32_t>&> {
				if constexpr (reverse)
					return EDG_rev<uint32_t>(&eget<ED_Graph<uint32_t>>(env));
				else
					return eget<ED_Graph<uint32_t>>(env);
			}(),
			EDG_view<uint32_t, reverse>(&eget<ED_Graph<uint32_t>>(env)),
			events, events_rev,
			left_mate,
			tmp_alignments_l, MAX_ALIGNS_M[pair_idx],
			[&] (bool r_c, uint32_t i_r, double val) {
				const auto& [eds_pos, read_pos] = (r_c ? pos_pairs_a_r_c : pos_pairs_b_r_c)[i_r];
				return hotspot<int_type>(eds_pos, read_pos, (int_type) std::max(0.0, PE_FRAGMENT_LENGTH - (double)val));
			},
			[&] (std::vector<hotspot<int_type>>&& hotspots_l, std::vector<hotspot<int_type>>&& hotspots_c) {
				return read_processor::align<int_type>
					( std::move(hotspots_l), read_b_r_c
					, std::move(hotspots_c), read_a_r_c
					, env, MAX_DIST[align_idx]
					, MAX_ALIGNS_C[align_idx], MAX_ALIGNS_T[align_idx]
					);
			},
			mk_event_pos);
		if (not tmp_alignments_r.empty()) break;
	}

	if (not tmp_alignments_r.empty())
	{ // we have found a mate pair
		std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> mates; // score, tmp_l, tmp_r
		for (size_t i = 0; i < tmp_alignments_r.size(); i++) {
			const auto src = left_mate[ tmp_alignments_r[i].reverse_compl ][ tmp_alignments_r[i].source ];
			assert(tmp_alignments_l[src].reverse_compl == tmp_alignments_r[i].reverse_compl);

			mates.emplace_back( tmp_alignments_l[src].get_dist() + tmp_alignments_r[i].get_dist(),
					src, i);
		}
		std::sort(mates.begin(), mates.end(), [](const auto& lhs, const auto& rhs) {
			return std::get<0>(lhs) < std::get<0>(rhs);
		});

		if (mates.size() > MAX_ALIGNS_O)
			mates.resize(MAX_ALIGNS_O);
		std::vector<temp_alignment<int_type>> tmp_l(mates.size()), tmp_r(mates.size());
		for (size_t i = 0; i < mates.size(); i++) {
			tmp_l[i] = tmp_alignments_l[std::get<1>(mates[i])];
			tmp_r[i] = tmp_alignments_r[std::get<2>(mates[i])];
		}
		return {
			read_processor::finalize_alignments(read_a, read_b, env, tmp_l, MAX_ALIGNS_O),
			read_processor::finalize_alignments(read_b_r_c, read_a_r_c, env, tmp_r, MAX_ALIGNS_O)
		};
	}
	else
	{ // we have not found mate pairs
		if constexpr (not reverse)
		{ // fallback
			if (not FALLBACK) // no fallback
				return {{{}, {}}};

			fallback_tries = 1;
			auto tmp_alignments = fallback<int_type>(
				env,
				events, events_rev,
				hotspots_a, hotspots_b,
				pos_pairs_a_r_c, pos_pairs_b_r_c,
				read_a, read_b,
				fallback_tries, pairing_time);
			std::sort(tmp_alignments.begin(), tmp_alignments.end(), [](const auto& lhs, const auto& rhs) {
				return lhs[0].get_dist() + lhs[1].get_dist() < rhs[0].get_dist() + rhs[1].get_dist();
			});
			std::vector< temp_alignment<int_type> > tmp_left, tmp_right;
			tmp_left.resize(std::min<size_t>(MAX_ALIGNS_O, tmp_alignments.size())), tmp_right.resize(tmp_left.size());
			for (size_t i = 0; i < tmp_left.size(); i++)
				std::tie(tmp_left[i], tmp_right[i]) = std::tie(tmp_alignments[i][0], tmp_alignments[i][1]);
			return {
				read_processor::finalize_alignments(read_a, read_b, env, tmp_left, MAX_ALIGNS_O),
				read_processor::finalize_alignments(read_b_r_c, read_a_r_c, env, tmp_right, MAX_ALIGNS_O)
			};
		}
		else
		{ // try matching right alignments to left query_positions
			return try_pair<int_type, not reverse>(
				env,
				read_a_r_c,
				read_b_r_c,
				read_a,
				read_b,
				pos_pairs_a_r_c,
				pos_pairs_b_r_c,
				pos_pairs_a,
				pos_pairs_b,
				fallback_tries,
				pairing_time);
		}
	}
}

// align left part of fragment to reference ("left" according to the reference)
// NOTE: align_r_c = true means that the entire fragment is the reverse-complement
template<typename int_type, typename Environ>
std::array< std::vector<alignment<int_type>>, 2 >
pair_read(
	const Environ& env,
	const fasta_read<int_type>& read_l,
	const fasta_read<int_type>& read_r,
	size_t& fallback_tries,
	double& pairing_time
) {
	const auto& mini = eget<gedmap_mini::minimizer_index>(env);

	const auto get_positions = [&](const fasta_read<int_type>& read) {
		return read_processor::get_positions<int_type>(
			read_processor::get_fragments<int_type>(read, FRAGMENT_COUNT[0], mini),
			mini);
	};

	const auto& read_l_r_c = read_l.get_rev_compl();
	const auto& read_r_r_c = read_r.get_rev_compl();

	std::vector< pe_event_type > events, events_rev;

	auto pos_pairs_l = get_positions(read_l);
	auto pos_pairs_r = get_positions(read_r);
	auto pos_pairs_r_r_c = get_positions(read_r_r_c);
	auto pos_pairs_l_r_c = get_positions(read_l_r_c);

	return try_pair<int_type, true>(
		env,
		read_l_r_c,
		read_r_r_c,
		read_l,
		read_r,
		pos_pairs_l_r_c,
		pos_pairs_r_r_c,
		pos_pairs_l,
		pos_pairs_r,
		fallback_tries,
		pairing_time);
}



// NOTE: we assume FR-reads, i.e. one read aligns to the forward strand and the
// other aligns to the reverse strand at a later position so that they point
// towards each other. We of course also consider the case that the fragment is
// the reverse-complement
template< class int_type >
void
map_pairs(
	gedmap_mini::minimizer_index & mini,
	const ED_Graph<uint32_t>& msq, // minimum shift query
	const std::string & EDS,
	const adjacency & adj,
	const pos_EDS_to_FA_type & p2FA,
	std::istream& fastq_l,
	std::istream& fastq_r,
	std::ofstream& o_s)
{
	environment<
		gedmap_mini::minimizer_index,
		std::string, // EDS
		adjacency,
		pos_EDS_to_FA_type,
		ED_Graph<uint32_t>
		> env(mini, EDS, adj, p2FA, msq);
	// set large buffers for fastq_{l,r}. Mostly useful for many threads
	std::vector<char> outbuffer(1<<23);
	o_s.rdbuf()->pubsetbuf(outbuffer.data(), outbuffer.size());

	gedmap_io::print_row("Fragment length ", PE_FRAGMENT_LENGTH);
	if(THREAD_COUNT) omp_set_num_threads(THREAD_COUNT);

	const auto num_reads = [&] {
		std::string tmp;
		size_t num_lines = 0;
		while (std::getline(fastq_l, tmp))
			num_lines++;
		fastq_l.clear();
		fastq_l.seekg(0, std::ios::beg);
		return num_lines / 4;
	} ();

	size_t num_failed = 0;
	size_t num_fallback = 0, num_fallback_tries = 0;
	double pairing_time = 0;

	size_t count = 0;
	#pragma omp parallel for schedule(dynamic,10)
	for (size_t read_i = 0; read_i < num_reads; read_i++) {
		fasta_read<int_type> read_l, read_r;
		#pragma omp critical
		{
			read_l = fasta_read<int_type>(fastq_l);
			read_r = fasta_read<int_type>(fastq_r);
		}

		size_t fallback_tries = 0;
		double pairing_time_tmp = 0;
		auto[alignments_l, alignments_r] = pair_read<int_type>(env, read_l, read_r, fallback_tries, pairing_time_tmp);

		#pragma omp critical
		{
			const auto cnt = read_processor::write_aligned_mates(
				alignments_l, alignments_r,
				read_l, read_r,
				o_s,
				p2FA);

			if (fallback_tries > 0) {
				num_fallback++;
				num_fallback_tries += fallback_tries - 1;
			}
			pairing_time += pairing_time_tmp;

			if (cnt == 0)
				num_failed++;
			count++;
			if (count % 1000 == 0)
				gedmap_io::flush_row("Searched read pairs", to_string(count));
		}
		
	}
	std::cerr << std::endl << num_failed << " could not pair" << std::endl;
	std::cerr << "num_fallback: " << num_fallback << std::endl;
	std::cerr << "num_fallback_tries: " << num_fallback_tries << std::endl;
	std::cerr << "fallback pairing time: " << pairing_time << "s" << std::endl;
	std::cerr << "paired_end_skip_queries: " << paired_end_skip_queries.load() << std::endl;
	std::cerr << "paired_end_skip: " << paired_end_skip.load() << std::endl;
}

} // namespace gedmap_align_min
