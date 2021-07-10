from dash import Dash
from dash.dependencies import Input, Output, State
from dash_html_components import Div, H1, Hr, Button, Label, H3, Table, Tr, Td, Pre
# from dash_marvinjs import DashMarvinJS
from dash_table import DataTable
from dash_core_components import Dropdown, Input as InputField, Markdown, Slider, Graph, Tabs, Tab
from pickle import load
from pandas import DataFrame
import pandas as pd
from CGRtools.files import RDFRead, RDFWrite
from CGRtools import ReactionContainer
import dash_dangerously_set_inner_html as dhtml
import pprint


external_stylesheets = [{'href': 'https://stackpath.bootstrapcdn.com/bootstrap/4.1.3/css/bootstrap.min.css',
                         'rel': 'stylesheet', 'crossorigin': 'anonymous',
                         'integrity': 'sha384-MCw98/SFnGE8fJT3GXwEOngsV7Zt27NXFoaoApmYm81iuXoPkFOJwJ8ERdknLPMO'}]
external_scripts = [{'src': 'https//cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.4/MathJax.js?config=TeX-MML-AM_CHTML'}]
app = Dash(__name__, external_stylesheets=external_stylesheets, external_scripts=external_scripts)
desc = '''
  __Instruction for collections remapping rules:__
* Choose the reaction center in table;
* Click on "reaction center" or "reaction" to identify the AAM;
* If the AAM is correct, click on "Its correct" and reaction center will be saved as correct;
* If the AAM is wrong, click on "Its wrong", add the correct and incorrect atom numbers of the products to the table and click "Save remapping rule". '''
text = Div([Markdown(desc)])

df = pd.read_csv('table_for_web.csv')
df = df[['reaction center', 'frequency of rc']]
with open('reaction_center_uspto_end.pickle', 'rb') as f:
    reaction_centers = load(f)

tabs = Tabs(id='tabs', value='good',
            children=[Tab(label='reaction center', value='reaction center', id='reaction center'),
                        Tab(label='reaction', value='reaction', id='reaction')])
tabs_out = Div(id='tabs-example-content')

table = Div([DataTable(id='table', columns=[{"name": i, "id": i} for i in df.columns], data=df.to_dict('records'),
            row_selectable="single",
            style_header=dict(backgroundColor="paleturquoise"), style_data=dict(backgroundColor="lavender", whiteSpace='normal'),
            style_table={'maxHeight': '500px', 'overflowY': 'scroll'},
           style_cell_conditional=[ {'if': {'column_id': 'reaction center'}, 'width': '70%'},
            {'if': {'column_id': 'frequency of rc'}, 'width': '20%'}])],
            className='col-md-5')

table_fix = Div([
    DataTable(id='fix_numbers',
            columns=[{ 'name': 'Num {}'.format(i), 'id': 'Num-{}'.format(i)} for i in ('false', 'true')],
            data=[{'Num-{}'.format(i): (0,) for i in ('false', 'true')}],
            editable=True, row_deletable=True),
    Button('Add row', id='add-rows-button', n_clicks=0),
    Div(id='fix_numbers_output')
])

question = Div([Label(["What do you think about this mapping?", Button('Its correct', id='button_true_aam', n_clicks=0),
                                                                Button('Its wrong', id='button_false_aam', n_clicks=0)]),
                Div(id='out_true_button'),
                Div(id='out_false_button')])
clicks = Div([
            Hr(),
            table_fix,
            Hr(),
            Div([Label([Button('Save remmaping rule', id='button_rule', n_clicks=0)]),
                 Div(id='out_rule')])])
tabs_all_clicks = Div([tabs, tabs_out, question, clicks], className='col-md-7')
up = Div([table, tabs_all_clicks], className='row col-md-12')
app.layout = Div([H1("Collection remmaping rules", style={'textAlign': 'center'}), up, Hr(), text])


@app.callback(Output('tabs-example-content', 'children'),
              [Input('tabs', 'value')],
              [State("table", "selected_rows")])
def output_reaction_rc(tab, selected_rows):
    if selected_rows and tab == 'reaction':
        with RDFRead('reactions_uspto_for_web_end.rdf') as f:
            for num, r in enumerate(f):
                if num == selected_rows[0]:
                    return dhtml.DangerouslySetInnerHTML(r.depict())

    elif selected_rows and tab == 'reaction center':
        return dhtml.DangerouslySetInnerHTML(reaction_centers[selected_rows[0]].depict())

reaction_centers_true, reaction_centers_false = [], []

@app.callback(Output('out_true_button', 'children'),
            [Input('button_true_aam', 'n_clicks')],
            [State("table", "selected_rows")])
def correct_reaction(button, selected_rows):
    if selected_rows:
        rc = reaction_centers[selected_rows[0]]
        # reaction_centers_true.append(rc)
        return f'True reaction center {str(rc)} was saved. All true rc: {button}'

@app.callback(Output('out_false_button', 'children'),
            [Input('button_false_aam', 'n_clicks')],
            [State("table", "selected_rows")])
def wrong_reaction(button, selected_rows):
    if selected_rows:
        rc = reaction_centers[selected_rows[0]]
        # reaction_centers_false.append(rc)
        return f'False reaction center {str(rc)} was saved. All false rc: {button}'

@app.callback(Output('fix_numbers_output', 'children'), [Input('fix_numbers', 'data')])
def display_output(rows):
    return Div([
        Div('Fixing numbers'),
        Pre(pprint.pformat(rows))])


@app.callback(Output('fix_numbers', 'data'), [Input('add-rows-button', 'n_clicks')],
             [State('fix_numbers', 'data')], [State('fix_numbers', 'columns')])
def add_row(n_clicks, rows, columns):
    if n_clicks > 0:
        rows.append({c['id']: '' for c in columns})
    return rows


def load_rule(reactions):
    rules = []
    for bad, good in reactions:
        if str(bad) != str(good):
            raise ValueError('bad and good reaction should be equal')

        cgr_good, cgr_bad = ~good, ~bad
        gc = cgr_good.augmented_substructure(cgr_good.center_atoms, deep=1)
        bc = cgr_bad.augmented_substructure(cgr_bad.center_atoms, deep=1)
        atoms = set(bc.atoms_numbers + gc.atoms_numbers)

        re_g, re_b, pr_g, pr_b = set(), set(), set(), set()
        for pr in good.reactants:
            re_g.update(pr)
        for pr in good.products:
            pr_g.update(pr)
        for pr in bad.products:
            pr_b.update(pr)
        for pr in bad.reactants:
            re_b.update(pr)
        atoms.update((re_b.difference(pr_b)).intersection(pr_g))

        strange_atoms = pr_b.difference(pr_g)
        atoms.update(strange_atoms)
        bad_query = (cgr_bad).substructure(atoms.intersection(cgr_bad), as_query=True)
        good_query = (cgr_good).substructure(atoms.intersection(cgr_good), as_query=True)

        fix = {}
        for mb, mg in zip(bad.products, good.products):
            fix.update({k: v for k, v in zip(mb, mg) if k != v and k in atoms})  # get fix map
        valid = set(fix).difference(strange_atoms)
        rules.append((bad_query, good_query, fix, valid))
    return rules

RULES = []

@app.callback(Output('out_rule', 'children'),
            [Input('button_rule', 'n_clicks')],
            [Input('fix_numbers', 'data')],
            [State("table", "selected_rows")])
def save_rule(button, rows, selected_rows):
    if selected_rows:
        if button > 0:
            with RDFRead('reactions_uspto_for_web_end.rdf') as f:
                for num, r in enumerate(f):
                    if num == selected_rows[0]:
                        break
            mapping = dict()
            for d in rows:
                nums = [eval(i) for i in d.values()]
                mapping[nums[0]] = nums[1]
            mols_p = [p.remap(mapping, copy=True) for p in r.products]
            reaction_remap = ReactionContainer(r.reactants, mols_p, r.reagents)

            rule = load_rule([(r, reaction_remap)])
            fix = rule[0][2]

            if rule not in RULES and len(rule) > 0:
                RULES.append(rule)
            else:
                return f'Rule with fixing numbers {fix} is already saved or error. All saved rules: {len(RULES)}'
            return f'Rule with fixing numbers {fix} saved. All saved rules: {len(RULES)}'

if __name__ == '__main__':
    app.run_server(host='0.0.0.0', debug=False)
