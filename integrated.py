#!/usr/bin/env python3
# -*- coding: utf-8 -*
from __future__ import unicode_literals

import cgi
import pandas as pd
import numpy as np
from scipy.interpolate import griddata
from bs4 import BeautifulSoup
import requests
import plotly.graph_objs as go



gettext = cgi.FieldStorage()

particle_class=gettext.getfirst("Particle", "Pin")
w_start=float(gettext.getfirst("w_min", "0.1").replace(",", "."))
w_finish=float(gettext.getfirst("w_max", "4.0").replace(",", "."))
q_start=float(gettext.getfirst("q2_min", "0.5").replace(",", "."))
q_finish=float(gettext.getfirst("q2_max", "0.5").replace(",", "."))
energy=float(gettext.getfirst("energy", "5.75").replace(",", "."))
interpolation_step=float(gettext.getfirst("step", "0.01"))


particle_form_text=''
if particle_class == 'Pin':
    particle_form_text="""<p> Reaction channel:
                        <select class="select" name="Particle" size="1">
                        <option value="Pin">gvp--->π⁺n</option>
                        <option value="Pi0P">gvp--->π⁰p</option>
                        </select>
                    </p>"""
else:
    particle_form_text = """<p> Reaction channel:
                            <select class="select" name="Particle" size="1">
                            <option value="Pi0P">gvp--->π⁰p</option>
                            <option value="Pin">gvp--->π⁺n</option>
                            </select>
                        </p>"""

# particle_class="Pin"
# w_start=float("1.3")
# w_finish=float("1.3")
# q_start=float("0.5")
# q_finish=float("5")
# energy=float("5.75")



x_axis_name=''

mp = 0.93827
df = pd.read_csv('final_table.csv', header=None, sep='\t',
                     names=['Channel', 'MID', 'Wmin', 'Wmax', 'Q2min', 'Q2max', 'Cos(theta)', 'sigma_t', 'd_sigma_t',
                            'sigma_l', 'd_sigma_l', 'sigma_tt', 'd_sigma_tt', 'sigma_lt', 'd_sigma_lt', 'eps'])
df['w_average'] = (df['Wmin'] + df['Wmax']) / 2
df['q2_average'] = (df['Q2min'] + df['Q2max']) / 2




# read file and use data for particle
if particle_class == 'Pin':
    partNum = '1212'
    ParticleSecret = 'PIN'
    ParticleBeauty = 'gvp->π⁺n'
    dataframe = df[
        (df.Channel == 8) | (df.Channel == 14) | (df.Channel == 41) | (df.Channel == 141)].copy()
    dataframes = [
        dataframe[(dataframe['w_average'] >= 1.1) & (dataframe['w_average'] <= 1.6) &
                  (dataframe['q2_average'] >= 0.2) & (dataframe['q2_average'] <= 0.7)].copy(),
        dataframe[(dataframe['w_average'] >= 1.1) & (dataframe['w_average'] <= 1.15) &
                  (dataframe['q2_average'] >= 2.115) & (dataframe['q2_average'] <= 4.155)].copy(),
        dataframe[(dataframe['w_average'] >= 1.15) & (dataframe['w_average'] <= 1.69) &
                  (dataframe['q2_average'] >= 1.72) & (dataframe['q2_average'] <= 4.16)].copy(),
        dataframe[(dataframe['w_average'] >= 1.605) & (dataframe['w_average'] <= 2.01) &
                  (dataframe['q2_average'] >= 1.8) & (dataframe['q2_average'] <= 4)].copy()
    ]

elif particle_class == 'Pi0P':
    PartNum = '1213'
    ParticleSecret = 'PI0P'
    ParticleBeauty = 'gvp->π⁰p'
    dataframe = df[(df.Channel == 9) | (df.Channel == 37) | (df.Channel == 170)].copy()
    dataframes = [
        dataframe[(dataframe['w_average'] >= 1) & (dataframe['w_average'] <= 1.8) &
                  (dataframe['q2_average'] >= 0.3) & (dataframe['q2_average'] <= 1.9)],
        dataframe[(dataframe['w_average'] >= 1) & (dataframe['w_average'] <= 1.4) &
                  (dataframe['q2_average'] >= 2.9) & (dataframe['q2_average'] <= 6.1)]
    ]

# initialize method
method=0
if w_start==w_finish:
    x_axis_name='Q2 (GeV2)'
    method=1
    w_steps=[float(w_start)]
    # q_steps=np.linspace(float(q_start), float(q_finish), 200)
    q_steps= np.arange(float(q_start), float(q_finish)+interpolation_step, interpolation_step).tolist()
else:
    x_axis_name='W (GeV)'
    method=2
    q_steps=[float(q_start)]
    # w_steps=np.linspace(float(w_start), float(w_finish), 200)
    w_steps = np.arange(float(w_start), float(w_finish)+interpolation_step, interpolation_step).tolist()

cos_steps=np.linspace(-1,1,200)
w_values, q_values, cos_values = np.meshgrid(w_steps, q_steps, cos_steps, indexing='ij')
w_values, q_values, cos_values = w_values.ravel(), q_values.ravel(), cos_values.ravel()

# make interpolation and drop nans
sigma_t=[]
d_sigma_t=[]
sigma_l=[]
d_sigma_l=[]
sigma_tt=[]
d_sigma_tt=[]
sigma_lt=[]
d_sigma_lt=[]

w_res_val=[]
q_res_val=[]
cos_res_val=[]

for current_df in dataframes:
    points=(current_df.w_average.ravel().flatten(),
            current_df.q2_average.ravel().flatten(),
            current_df['Cos(theta)'].ravel().flatten())

    values=(current_df.sigma_t.ravel().flatten(),
           current_df.d_sigma_t.ravel().flatten(),
           current_df.sigma_l.ravel().flatten(),
           current_df.d_sigma_l.ravel().flatten(),
           current_df.sigma_tt.ravel().flatten(),
           current_df.d_sigma_tt.ravel().flatten(),
           current_df.sigma_lt.ravel().flatten(),
           current_df.d_sigma_lt.ravel().flatten())

    values=np.stack((values), axis=-1)

    current_interpolation = griddata(points, values,
                                     (w_values, q_values, cos_values),
                                     method='linear', rescale=True)

    not_nans = ~np.isnan(current_interpolation[:,0])

    sigma_t+=current_interpolation[:,0][not_nans].flatten().tolist()
    d_sigma_t+=current_interpolation[:,1][not_nans].flatten().tolist()
    sigma_l+=current_interpolation[:,2][not_nans].flatten().tolist()
    d_sigma_l+=current_interpolation[:,3][not_nans].flatten().tolist()
    sigma_tt+=current_interpolation[:,4][not_nans].flatten().tolist()
    d_sigma_tt+=current_interpolation[:,5][not_nans].flatten().tolist()
    sigma_lt+=current_interpolation[:,6][not_nans].flatten().tolist()
    d_sigma_lt+=current_interpolation[:,7][not_nans].flatten().tolist()

    w_res_val+=w_values[not_nans].flatten().tolist()
    q_res_val+=q_values[not_nans].flatten().tolist()
    cos_res_val+=cos_values[not_nans].flatten().tolist()

# sort values and calculate eps
tmp_df=pd.DataFrame({'w':w_res_val,
                     'q':q_res_val,
                     'cos':cos_res_val,
                     'sigma_t': sigma_t,
                     'd_sigma_t':d_sigma_t,
                     'sigma_l': sigma_l,
                     'd_sigma_l':d_sigma_l,
                     'sigma_tt': sigma_tt,
                     'd_sigma_tt':d_sigma_t,
                     'sigma_lt': sigma_lt,
                     'd_sigma_lt':d_sigma_lt})

tmp_df.sort_values(by=['w','q','cos'],inplace=True)
tmp_df.reset_index(drop=True,inplace=True)
tmp_df['nu']=(tmp_df['w'] ** 2 + tmp_df['q'] - mp * mp) / (2 * mp)
tmp_df['eps']=1/(1 + 2*(tmp_df['nu']**2 + tmp_df['q'])/(4*(energy - tmp_df['nu'])*energy-tmp_df['q']))
tmp_df['sigma_u']=tmp_df['sigma_t']+tmp_df['sigma_l']*tmp_df['eps']
tmp_df['sigma']=tmp_df['sigma_u']*2*np.pi

# calculate integral cross section
w_graph=[]
q_graph=[]
sigma_integrated=[]
d_sigma_integrated=[]

if method==2:
    abscissa='x'

    for el in tmp_df['w'].unique():
        t=tmp_df[tmp_df['w']==el].copy()
        t.reset_index(inplace=True, drop=True)
        integral_sigma=0
        d_integral_sigma=0
        for cos_idx in range(0,len(t)-1):
            integral_sigma += (t.loc[cos_idx+1, 'cos'] - t.loc[cos_idx, 'cos']) * \
                                ((t.loc[cos_idx+1, 'sigma'] + t.loc[cos_idx, 'sigma'])/2)
            d_integral_sigma += (t.loc[cos_idx+1, 'sigma']**2 + t.loc[cos_idx, 'sigma']**2) * \
                                 ((t.loc[cos_idx+1, 'cos'] - t.loc[cos_idx, 'cos'])/2)**2
        w_graph.append(el)
        q_graph.append(q_start)
        sigma_integrated.append(integral_sigma)
        d_sigma_integrated.append(d_integral_sigma**0.5)
        

if method==1:
    abscissa='Q2'
    for el in tmp_df['q'].unique():
        t=tmp_df[tmp_df['q']==el].copy()
        t.reset_index(inplace=True, drop=True)
        integral_sigma=0
        d_integral_sigma=0
        for cos_idx in range(0,len(t)-1):
            integral_sigma += (t.loc[cos_idx+1, 'cos'] - t.loc[cos_idx, 'cos']) * \
                                ((t.loc[cos_idx+1, 'sigma'] + t.loc[cos_idx, 'sigma'])/2)
            d_integral_sigma += (t.loc[cos_idx+1, 'sigma']**2 + t.loc[cos_idx, 'sigma']**2) * \
                                 ((t.loc[cos_idx+1, 'cos'] - t.loc[cos_idx, 'cos'])/2)**2
        w_graph.append(w_start)
        q_graph.append(el)
        sigma_integrated.append(integral_sigma)
        d_sigma_integrated.append(d_integral_sigma**0.5)


set_channel =  [255,1] if particle_class=='Pin' else [255,2]
set_channel_name = 'π+n' if particle_class=='Pin' else 'π0p'



url='https://clas.sinp.msu.ru/strfun/'
values = {'quantity':'sigma',
          'RLT-src':'',
          'channel':set_channel,
          'q2-first':q_start,
          'q2-step':interpolation_step,
          'q2-last':q_finish,
          'abscissa':abscissa,
          'x-first':w_start,
          'x-step':interpolation_step,
          'x-last':w_finish,
          'ebeam':5.75,
          'L':12.8e10,
          'DeltaW':interpolation_step,
          'DeltaQ2':interpolation_step,
          'dataset':0,
          'res-view':'html',
          'linetype':'lines',
          'submit':'Calculate'}

response = requests.get(url, params=values)
if response.ok:
    soup = BeautifulSoup(response.text, 'html.parser')
    table = soup.find('table', {'class': 'sortable'})
    df_vitaly = pd.read_html(str(table))[0]
else:
    pass
#     print(f'Request failed with status code {response.status_code}')




def graph_maker(x_array=[], y_array=[], d_y_array=[],
                         x_exp_data=[], y_exp_data=[], dy_exp_data=[],
                         x_exp_data_inclusive=[], y_exp_data_inclusive=[], dy_exp_data_inclusive=[],
                         layout_title='Integral cross section (mcbn)', x_label=x_axis_name):

    trace_interp = go.Scatter(
        x=x_array,
        y=y_array,
        error_y=dict(
            type='data',
            array=d_y_array,
#             color='rgba(100, 100, 255, 0.6)',
            thickness=1.5,
            width=3),
        name='interpolation Almaz ' + str('π+n' if particle_class=='Pin' else 'π0p'),
        marker_size=1)

    trace_interp_vitaly = go.Scatter(
        x=x_exp_data,
        y=y_exp_data,
        error_y=dict(
            type='data',
            array=dy_exp_data,
#             color='rgba(100, 20, 255, 0.6)',
            thickness=1.5,
            width=3),
        name='interpolation Vitaly ' + str('π+n' if particle_class=='Pin' else 'π0p'),
        marker_size=1)


    trace_interp_inclusive = go.Scatter(
#         mode='markers',
        x=x_exp_data_inclusive,
        y=y_exp_data_inclusive,
        name='interpolation Vitaly inclusive',
        marker=dict(color='rgba(100, 100, 100, 1)',
                    symbol='square'),
        error_y=dict(
            type='data',
            array=dy_exp_data_inclusive,
#             color='rgba(100, 100, 100, 1)',
            thickness=1.5,
            width=3),
        marker_size=10)

    data = [trace_interp, trace_interp_vitaly, trace_interp_inclusive]

    fig = go.Figure(data=data)
    fig.layout.height = 700
    fig.layout.width = 1000
    fig.layout.title = layout_title

    fig.layout.yaxis = dict(
        showgrid=True,
        zeroline=True,
        showline=True,
        gridcolor='#bdbdbd',
        gridwidth=1,
        zerolinecolor='black',
        zerolinewidth=0.5,
        linewidth=0.5,
        title=layout_title,
        titlefont=dict(
            family='Arial, sans-serif',
            size=18,
            color='black'
        ))
    fig.layout.xaxis = dict(
        showgrid=True,
        zeroline=True,
        showline=True,
        gridcolor='#bdbdbd',
        gridwidth=1,
        zerolinecolor='black',
        zerolinewidth=0.5,
        linewidth=0.2,
        title=x_label,
        titlefont=dict(
            family='Arial, sans-serif',
            size=18,
            color='black'
        ))

    return fig




def make_part_graph(x_array=[], y_array=[], d_y_array=[],
                 y_exp_data=[], dy_exp_data=[],
                layout_title='Part of the inclusive cross section', x_label=x_axis_name):
    
    trace_interp = go.Scatter(
        x=x_array,
        y=y_array,
        error_y=dict(
            type='data',
            array=d_y_array,
#             color='rgba(100, 100, 255, 0.6)',
            thickness=1.5,
            width=3),
        name='interpolation Almaz ' + str('π+n' if particle_class=='Pin' else 'π0p'),
        marker_size=1)

    trace_interp_vitaly = go.Scatter(
        x=x_array,
        y=y_exp_data,
        error_y=dict(
            type='data',
            array=dy_exp_data,
#             color='rgba(100, 20, 255, 0.6)',
            thickness=1.5,
            width=3),
        name='interpolation Vitaly ' + str('π+n' if particle_class=='Pin' else 'π0p'),
        marker_size=1)

    data = [trace_interp, trace_interp_vitaly]

    fig = go.Figure(data=data)
    fig.layout.height = 700
    fig.layout.width = 1000
    fig.layout.title = layout_title

    fig.layout.yaxis = dict(
        showgrid=True,
        zeroline=True,
        showline=True,
        gridcolor='#bdbdbd',
        gridwidth=1,
        zerolinecolor='black',
        zerolinewidth=0.5,
        linewidth=0.5,
        title=layout_title,
        titlefont=dict(
            family='Arial, sans-serif',
            size=18,
            color='black'
        ))
    fig.layout.xaxis = dict(
        showgrid=True,
        zeroline=True,
        showline=True,
        gridcolor='#bdbdbd',
        gridwidth=1,
        zerolinecolor='black',
        zerolinewidth=0.5,
        linewidth=0.2,
        title=x_label,
        titlefont=dict(
            family='Arial, sans-serif',
            size=18,
            color='black'
        ))

    return fig



if method==1:
    fig=graph_maker(x_array=q_graph,
                             y_array=sigma_integrated,
                             d_y_array=d_sigma_integrated,
                         x_exp_data=df_vitaly[df_vitaly['Channel']==set_channel_name]['Q2, GeV2'],
                         y_exp_data=df_vitaly[df_vitaly['Channel']==set_channel_name]['σ, μb'],
                         dy_exp_data=df_vitaly[df_vitaly['Channel']==set_channel_name]['Δσ, μb'],

                         x_exp_data_inclusive=df_vitaly[df_vitaly['Channel']=='inclusive']['Q2, GeV2'],
                         y_exp_data_inclusive=df_vitaly[df_vitaly['Channel']=='inclusive']['σ, μb'],
                         dy_exp_data_inclusive=df_vitaly[df_vitaly['Channel']=='inclusive']['Δσ, μb'],
                         layout_title='Integral cross section (mcbn)', x_label=x_axis_name)
    

else:
    fig=graph_maker(x_array=w_graph,
                             y_array=sigma_integrated,
                             d_y_array=d_sigma_integrated,
                         x_exp_data=df_vitaly[df_vitaly['Channel']==set_channel_name]['W, GeV'],
                         y_exp_data=df_vitaly[df_vitaly['Channel']==set_channel_name]['σ, μb'],
                         dy_exp_data=df_vitaly[df_vitaly['Channel']==set_channel_name]['Δσ, μb'],

                         x_exp_data_inclusive=df_vitaly[df_vitaly['Channel']=='inclusive']['W, GeV'],
                         y_exp_data_inclusive=df_vitaly[df_vitaly['Channel']=='inclusive']['σ, μb'],
                         dy_exp_data_inclusive=df_vitaly[df_vitaly['Channel']=='inclusive']['Δσ, μb'],
                         layout_title='Integral cross section (mcbn)', x_label=x_axis_name)

    
    
    
if method==1:
    min_val_second_interpolation=max(min(q_graph),
        df_vitaly[df_vitaly['Channel']=='inclusive']['Q2, GeV2'].min(),
        df_vitaly[df_vitaly['Channel']==set_channel_name]['Q2, GeV2'].min())
    max_val_second_interpolation=min(max(q_graph),
        df_vitaly[df_vitaly['Channel']=='inclusive']['Q2, GeV2'].max(),
        df_vitaly[df_vitaly['Channel']==set_channel_name]['Q2, GeV2'].max())


if method==2:
    min_val_second_interpolation=max(min(w_graph),
        df_vitaly[df_vitaly['Channel']=='inclusive']['W, GeV'].min(),
        df_vitaly[df_vitaly['Channel']==set_channel_name]['W, GeV'].min())
    max_val_second_interpolation=min(max(w_graph),
        df_vitaly[df_vitaly['Channel']=='inclusive']['W, GeV'].max(),
        df_vitaly[df_vitaly['Channel']==set_channel_name]['W, GeV'].max())


common_x_axis=np.arange(round(min_val_second_interpolation,2), 
                        round(max_val_second_interpolation,2),
                        0.001)

if method==1:
    sigma_integrated_grid = griddata(np.array(q_graph), 
                                      np.array(sigma_integrated),
                                      np.array(common_x_axis),
                                      method='linear', rescale=True)
    
    sigma_integrated_grid_vitaly = griddata(
                                      np.array(df_vitaly[df_vitaly['Channel']==set_channel_name]['Q2, GeV2']),
                                      np.array(df_vitaly[df_vitaly['Channel']==set_channel_name]['σ, μb']),
                                      np.array(common_x_axis),
                                      method='linear', rescale=True)

    sigma_integrated_grid_inclusive = griddata(
                                      np.array(df_vitaly[df_vitaly['Channel']=='inclusive']['Q2, GeV2']),
                                      np.array(df_vitaly[df_vitaly['Channel']=='inclusive']['σ, μb']),
                                      np.array(common_x_axis),
                                      method='linear', rescale=True)
    
    
if method==2:
    sigma_integrated_grid = griddata(np.array(w_graph), 
                                      np.array(sigma_integrated),
                                      np.array(common_x_axis),
                                      method='linear', rescale=True)
    
    sigma_integrated_grid_vitaly = griddata(
                                      np.array(df_vitaly[df_vitaly['Channel']==set_channel_name]['W, GeV']),
                                      np.array(df_vitaly[df_vitaly['Channel']==set_channel_name]['σ, μb']),
                                      np.array(common_x_axis),
                                      method='linear', rescale=True)

    sigma_integrated_grid_inclusive = griddata(
                                      np.array(df_vitaly[df_vitaly['Channel']=='inclusive']['W, GeV']),
                                      np.array(df_vitaly[df_vitaly['Channel']=='inclusive']['σ, μb']),
                                      np.array(common_x_axis),
                                      method='linear', rescale=True)

    
result_df=pd.DataFrame({'x_axis_values':common_x_axis,
              'sigma_integrated_grid' : sigma_integrated_grid,
              'sigma_integrated_grid_vitaly' : sigma_integrated_grid_vitaly,
              'sigma_integrated_grid_inclusive' : sigma_integrated_grid_inclusive })

result_df.dropna(inplace=True)
result_df.reset_index(inplace=True, drop=True)

result_df['frac_sigma']=result_df['sigma_integrated_grid']/result_df['sigma_integrated_grid_inclusive']
result_df['frac_sigma_vitaly']=result_df['sigma_integrated_grid_vitaly']/result_df['sigma_integrated_grid_inclusive']

fig_part = make_part_graph(x_array=result_df['x_axis_values'],
                           y_array=result_df['frac_sigma'],
                           d_y_array=[],
                           y_exp_data=result_df['frac_sigma_vitaly'],
                           dy_exp_data=[],
                           layout_title='Part of the inclusive cross section',
                           x_label=x_axis_name)



print("Content-type: text/html\n")
print("""<!DOCTYPE HTML>
<html>
<head>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style type="text/css">
        A {
            text-decoration: none;
            color: red;
        }

        * {
            margin: 0;
        }

        .textBox {
            width: 1440px;
            height: 80px;
            margin: auto;
        }

        .imagesBox {
            width: 1440px;
            height: 900px;
            margin: auto;
        }

        .textBox2 {
            width: 1440px;
            height: 50px;
            margin: auto;
        }

        .tableBox1 {
            margin: auto;
            width: 1440px;
            height: 350px;
        }

        .checkbox_msg {
            color: blue;
        }

        td {
            text-align: center;
        }

        .first_box {
            background-color: rgba(200, 200, 200, 0.6);
            width: 1070px;
            height: 570px;
            margin: auto;
            border-radius: 10px;
            margin-bottom: 30px;
        }


        .box_in_the_box {
            background-color: rgba(100, 100, 100, 0.6);
            position: absolute;
            width: 300px;
            height: 520px;
            margin-top: 25px;
            margin-left: 25px;
            border-radius: 4px;
        }

        .second_box_in_the_box {
            position: absolute;
            border-radius: 4px;
            margin-left: 350px;
            margin-top: 25px;
            width: 700px;
            height: 520px;
            background-color: rgba(100, 100, 100, 0.6);
        }

        .left_box {
            position: absolute;
            border-radius: 4px;
            margin-left: 20px;
            margin-top: 10px;
            width: 330px;
            height: 430px;
            padding: 10px;
            background-color: rgba(255, 255, 255, 0.6);

        }


        .right_box {
            position: absolute;
            border-radius: 4px;
            margin-left: 400px;
            margin-top: 10px;
            width: 280px;
            height: 450px;
            background-color: rgba(255, 255, 255, 0.6);

        }

        .input_small {
            width: 80px;
        }


        .left_sub_box {

            position: absolute;
            width: 130px;
            height: 280px;


        }

        .right_sub_box {
            margin-top: -40px;
            position: absolute;
            width: 130px;
            margin-left: 130px;
            height: 280px;
        }

        .select_hidden {
            visibility: hidden;
        }""")

print(""""
    </style>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width">
    <script type="text/javascript"
            src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
    </script>
    <title>CLAS graph</title>
</head>
<body>
<center>

    <br>
    <br>
    <br>
	<a href="https://clas.sinp.msu.ru/cgi-bin/almaz/instruction">Available data areas</a>
	<br>

    <br><br>


    <form method="GET" action="https://clas.sinp.msu.ru/cgi-bin/almaz/integrated/integrated.py">

        <div class="first_box">
            <div class="box_in_the_box">
                <br><br><br><br><br><br><br><br><br><br>
            </div>
            
            <div class="second_box_in_the_box">
                <br>
                Integrated cross section
                <div class="left_box">
                    <br>
                    <br>
                    Для получение интегрального сечения как функции W введите различные
                    w_min и w_max и одинаковое значение q_min=q_max и энергию реакции
                    <br>
                    (пример: w_min=0.4,  w_max=4,  q_min=0.5, q_max=0.5, E_beam=5.75)
                    <br>
                    ___________________________
                    <br><br>
                    Для получение интегрального сечения как функции Q2 введите различные
                    q_min и q_max и одинаковое значение w_min=w_max и энергию реакции
                    <br>
                    (пример: w_min=1.4,  w_max=1.4,  q_min=0.5, q_max=4, E_beam=5.75)
                    <br>
                    ___________________________
                    <br><br>
                </div>

                <div class="right_box">
                    <br><br><br><br>
                    <p>W min (GeV):&nbsp;&nbsp;&nbsp;&nbsp; <input type="text" class="input_small"  name="w_min" value="{}"
                                                                                placeholder="W min (GeV)"
                                                                                ></p>
                    <br>
                    <p>W max (GeV):&nbsp;&nbsp;&nbsp;&nbsp; <input type="text" class="input_small"  name="w_max" value="{}"
                                                                  placeholder="W max (GeV)"
                                                                ></p>
                    <br>
                    <p>Q2 min (GeV2):&nbsp;&nbsp;&nbsp;<input type="text" class="input_small"  name="q2_min" value="{}"
                                                                                   placeholder="Q2 min (GeV2)"
                                                                                   >
                    </p>
                    <br>
                    <p>Q2 max (GeV2):&nbsp;&nbsp;<input type="text" class="input_small" name="q2_max" value="{}"
                                                       placeholder="Q2 max (GeV2)"
                                                             ></p>

                    <br><br><br><p>
                        &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Beam energy:
                        <input class="input_small" type="text" name="eBeam" placeholder="MeV" value="{}">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
                    </p>
                    <br>

                    <p> Interpolation step:
                        <input class="input_small" type="text" name="grid_step_user" value="{}"
                               placeholder="grid step"> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
                    </p>
                    <br><br><br> 
                    
                    {}
                    
                    <br>
                    <p><input class="button" class="submitbutton" type="submit" value="Run"></p>

                </div>
            </div>
           nasrtdinov.ag17@physics.msu.ru
        </div>
    </form>""".format(w_start,w_finish,q_start,q_finish,energy,interpolation_step,particle_form_text))


particle_class=gettext.getfirst("Particle", "Pin")
w_start=float(gettext.getfirst("w_min", "0.1").replace(",", "."))
w_finish=float(gettext.getfirst("w_max", "4.0").replace(",", "."))
q_start=float(gettext.getfirst("q2_min", "0.5").replace(",", "."))
q_finish=float(gettext.getfirst("q2_max", "0.5").replace(",", "."))
energy=float(gettext.getfirst("energy", "5.75").replace(",", "."))
interpolation_step=float(gettext.getfirst("step", "0.01"))




print("""
        <br><br>
        Current values:
        <br><br>
        particle:  {}    &nbsp;&nbsp;&nbsp;&nbsp;   beam_energy: {}   &nbsp;&nbsp;&nbsp;&nbsp; interpolation step : {}
        <br><br>
        w_min = {} &nbsp;&nbsp;&nbsp;   w_max={}  &nbsp;&nbsp;&nbsp;
        q2_min={}   &nbsp;&nbsp;&nbsp;     q2_max={}""".format(particle_class,energy,interpolation_step,w_start,w_finish,q_start,q_finish))


print("{}".format(fig.to_html(full_html=False)))
print("<br><br><br><br>")
print("{}".format(fig_part.to_html(full_html=False)))
print("""</center>
        </body>
         </html>""")