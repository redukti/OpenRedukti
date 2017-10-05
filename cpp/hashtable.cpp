/**
 * DO NOT REMOVE COPYRIGHT NOTICES OR THIS HEADER.
 *
 * Contributor(s):
 *
 * The Original Software is OpenRedukti (https://github.com/redukti/OpenRedukti).
 * The Initial Developer of the Original Software is REDUKTI LIMITED (http://redukti.com).
 * Authors: Dibyendu Majumdar
 *
 * Copyright 2017 REDUKTI LIMITED. All Rights Reserved.
 *
 * The contents of this file are subject to the the GNU General Public License
 * Version 3 (https://www.gnu.org/licenses/gpl.txt).
 */

/**
 * Portions Copyright (C) 1994-2016 Lua.org, PUC-Rio.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 ******************************************************************************/
#include <hashtable.h>
#include <limits>

namespace redukti
{

template <> struct NilHelper<int> {
	static void set_nil(int *v) { *v = std::numeric_limits<int>::min(); }
	static bool is_nil(const int *v) { return *v == std::numeric_limits<int>::min(); }
};

template <> struct NilHelper<std::string> {
	static void set_nil(std::string *v) { v->clear(); }
	static bool is_nil(const std::string *v) { return v->length() == 0; }
};

int test_hashtable()
{
	typedef Table<int, int> TableIntInt;
	TableIntInt t;

	auto valp = t.newkey(5);
	valp->setval(66);
	if (t.get(5)->val() != 66)
		return 1;

	valp = t.newkey(3);
	valp->setval(67);
	if (t.get(5)->val() != 66 || t.get(3)->val() != 67)
		return 1;
	// t.dump();

	valp = t.newkey(7);
	valp->setval(68);
	if (t.get(5)->val() != 66 || t.get(3)->val() != 67 || t.get(7)->val() != 68)
		return 1;
	// t.dump();

	valp = t.newkey(2);
	valp->setval(69);
	if (t.get(5)->val() != 66 || t.get(3)->val() != 67 || t.get(7)->val() != 68 || t.get(2)->val() != 69)
		return 1;
	// t.dump();

	valp = t.newkey(11);
	valp->setval(70);
	if (t.get(5)->val() != 66 || t.get(3)->val() != 67 || t.get(7)->val() != 68 || t.get(2)->val() != 69 ||
	    t.get(11)->val() != 70)
		return 1;
	// t.dump();

	valp = t.newkey(19);
	valp->setval(71);
	if (t.get(5)->val() != 66 || t.get(3)->val() != 67 || t.get(7)->val() != 68 || t.get(2)->val() != 69 ||
	    t.get(11)->val() != 70 || t.get(19)->val() != 71)
		return 1;
	// t.dump();

	valp = t.newkey(4);
	valp->setval(72);
	if (t.get(5)->val() != 66 || t.get(3)->val() != 67 || t.get(7)->val() != 68 || t.get(2)->val() != 69 ||
	    t.get(11)->val() != 70 || t.get(19)->val() != 71 || t.get(4)->val() != 72)
		return 1;
	// t.dump();

	valp = t.newkey(14);
	valp->setval(73);
	if (t.get(5)->val() != 66 || t.get(3)->val() != 67 || t.get(7)->val() != 68 || t.get(2)->val() != 69 ||
	    t.get(11)->val() != 70 || t.get(19)->val() != 71 || t.get(4)->val() != 72 || t.get(14)->val() != 73)
		return 1;
	// t.dump();

	valp = t.newkey(-100);
	valp->setval(100);
	if (t.get(5)->val() != 66 || t.get(3)->val() != 67 || t.get(7)->val() != 68 || t.get(2)->val() != 69 ||
	    t.get(11)->val() != 70 || t.get(19)->val() != 71 || t.get(4)->val() != 72 || t.get(14)->val() != 73 ||
	    t.get(-100)->val() != 100)
		return 1;
	// t.dump();

	int count = 0;
	for (auto ptr = std::begin(t); ptr != std::end(t); ptr++) {
		// printf("%d %d\n", ptr->key(), ptr->value());
		count++;
	}
	if (count != 9) {
		return 1;
	}
	// printf("sizeof Node<int,int> = %d\n", (int)(sizeof(TableIntInt::node_type)));

	typedef Table<std::string, std::string> TableStrStr;
	TableStrStr t2;
	auto val = t2.newkey(std::string("hello"));
	val->setval("world!");
	val = t2.newkey(std::string("dibyendu"));
	val->setval("majumdar");
	if (t2.get("hello")->val() != std::string("world!") || t2.get("dibyendu")->val() != std::string("majumdar")) {
		return 1;
	}
	// t2.dump();

	// printf("sizeof Node<string,string> = %d, sizeof std::string() = %d\n",
	//       (int)(sizeof(TableStrStr::node_type)), (int)(sizeof(std::string)));
	printf("Hashtable Tests OK\n");
	return 0;
}
}