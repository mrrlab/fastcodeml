
#ifndef NEWICKGRAMMAR_H
#define NEWICKGRAMMAR_H

#include <boost/spirit/include/classic_core.hpp>
#include <boost/spirit/include/classic_parse_tree.hpp>

using namespace boost::spirit::classic;
////////////////////////////////////////////////////////////////////////////
//
//  Grammar for the Newick formatted tree files
//
////////////////////////////////////////////////////////////////////////////
struct NewickGrammar : public grammar<NewickGrammar>
{
    static const int treeID      = 10;
    static const int nodelistID  = 20;
    static const int subtreeID   = 30;
    static const int fulllabelID = 40;
    static const int branchlenID = 50;
    static const int cblenID     = 60;
    static const int labelID     = 70;
    static const int typeID      = 80;

    template <typename ScannerT>
    struct definition
    {
        definition(NewickGrammar const& /*self*/)
        {
            tree = nodelist >> !full_label >> !colon_plus_len >> no_node_d[*space_p >> ch_p(';')];

            nodelist = no_node_d[ch_p('(')] >> subtree % (no_node_d[ch_p(',') >> *space_p]) >> no_node_d[ch_p(')')];

            subtree = nodelist >> !full_label >> !colon_plus_len >> no_node_d[*space_p] | full_label >> !colon_plus_len >> no_node_d[*space_p];

            colon_plus_len = no_node_d[ch_p(':') >> *space_p] >> branch_length;

            branch_length = real_p;

            full_label = no_node_d[ch_p('#')] >> type | label >> !(no_node_d[ch_p('#')] >> type);

            label = leaf_node_d[alpha_p >> *(alnum_p|'.'|'_')];

            type = leaf_node_d[+alnum_p];
        }

        rule<ScannerT, parser_context<>, parser_tag<treeID> >      tree;
        rule<ScannerT, parser_context<>, parser_tag<nodelistID> >  nodelist;
        rule<ScannerT, parser_context<>, parser_tag<subtreeID> >   subtree;
        rule<ScannerT, parser_context<>, parser_tag<fulllabelID> > full_label;
        rule<ScannerT, parser_context<>, parser_tag<branchlenID> > branch_length;
        rule<ScannerT, parser_context<>, parser_tag<cblenID> >     colon_plus_len;
        rule<ScannerT, parser_context<>, parser_tag<labelID> >     label;
        rule<ScannerT, parser_context<>, parser_tag<typeID> >      type;

        rule<ScannerT, parser_context<>, parser_tag<treeID> >   const&
        start() const { return tree; }
    };
};

#endif


