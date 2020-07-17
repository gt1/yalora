/*
    yalora
    Copyright (C) 2018 German Tischler-Höhle

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#if ! defined(COORD_HPP)
#define COORD_HPP

#include <libmaus2/types/types.hpp>
#include <libmaus2/math/IntegerInterval.hpp>
#include <istream>
#include <sstream>

struct Coord
{
	uint64_t valid;
	uint64_t seq;
	uint64_t rc;
	uint64_t left;
	uint64_t length;

	static bool expect(std::istream & in, std::string const & ex)
	{
		uint64_t i = 0;

		while ( i < ex.size() )
		{
			int const c = in.get();

			if ( !in || c == std::istream::traits_type::eof() )
				return false;

			if ( c != ex[i] )
				return false;

			++i;
		}

		return i;
	}

	static bool getNumber(std::istream & in, uint64_t & v)
	{
		v = 0;
		uint64_t c = 0;

		while ( true )
		{
			if ( ! in )
			{
				if ( c )
					return true;
				else
					return false;
			}

			int const c = in.peek();

			if ( ! in || c == std::istream::traits_type::eof() || !isdigit(c) )
			{
				if ( c )
					return true;
				else
					return false;
			}

			v *= 10;
			v += c-'0';
			in.get();
		}
	}


	Coord()
	{

	}

	bool parse(std::string const & s)
	{
		std::istringstream istr(s);
		return parse(istr);
	}

	bool parse(std::istream & in)
	{
		bool ok = true;
		ok = ok && expect(in,"Coordinates(valid=");
		ok = ok && getNumber(in,valid);
		ok = ok && expect(in,",seq=");
		ok = ok && getNumber(in,seq);
		ok = ok && expect(in,",rc=");
		ok = ok && getNumber(in,rc);
		ok = ok && expect(in,",left=");
		ok = ok && getNumber(in,left);
		ok = ok && expect(in,",length=");
		ok = ok && getNumber(in,length);
		ok = ok && expect(in,")");
		return ok;
	}

	Coord merge(Coord const & C) const
	{
		Coord R;

		R.valid = valid;
		R.seq = seq;
		R.rc = rc;

		R.left = std::min(left,C.left);

		uint64_t const right = std::max(left+length,C.left+C.length);

		R.length = right - R.left;

		return R;
	}

	bool overlap(Coord const & C) const
	{
		if ( seq != C.seq )
			return false;
		if ( rc != C.rc )
			return false;

		libmaus2::math::IntegerInterval<int64_t> IA(
			left,
			left+length-1
		);
		libmaus2::math::IntegerInterval<int64_t> IC(
			C.left,
			C.left+C.length-1
		);

		return !IA.intersection(IC).isEmpty();
	}
};

std::ostream & operator<<(std::ostream & out, Coord const & C);
#endif
