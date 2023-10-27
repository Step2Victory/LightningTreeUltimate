import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd
import math, io, yaml
import asyncio
# import time, multiprocessing as mps


from enum import Enum
from plotly.subplots import make_subplots
from dash import Dash, dcc, html, Output, Input, State, ctx
from PIL import Image
import os
import subprocess
import flask

server = flask.Flask(__name__)

h_range = [4000, 12000]

path_to_project = '/home/' +os.environ['USER'] + '/LightningTreeUltimate/'

class SubprocessController:

    def __init__(self):
        self.task = None
        self.process = None
        self.event_loop = None

    def start(self, callback):
        self.stop() 
        self.task = asyncio.run(self.loop(callback))
        
    async def loop(self, callback):
        await self.internal_start()
        while await self.read():
            callback()
            self.write()

    async def internal_start(self):
        if self.process is None:
            print("Собираем бинарь плюсов")
            subprocess.run(["cmake","-DCMAKE_BUILD_TYPE=Release", "."])
            subprocess.run(["make"])
            print("Подпроцесс моделирования молнии запускается")
            self.process = await asyncio.create_subprocess_exec(path_to_project + "/main",
                                                stdout = asyncio.subprocess.PIPE,
                                                stdin = asyncio.subprocess.PIPE,
                                                stderr = asyncio.subprocess.STDOUT)
            print("Подпроцесс моделирования молнии запущен")
        else: print("Процесс уже запущен")
    
    def stop(self):
        if self.task is not None:
            self.task.cancel()
            self.process.terminate()
            self.process = None

    async def read(self):
        if self.process is not None:
            
            answer = await self.process.stdout.readline()
            print(answer)
            code, iter_number, charges_number = map(int, answer.split()[:-1])
            if code == 0:
                self.stop()
                return False

            time = float(answer[-1])
            print("Номер итерации: ", iter_number)
            print("Количество зарядов в графе", charges_number)
            print("Время от начала процесса", time)
            return True
        return False
    
    def write(self):
        if self.process is not None:
            self.process.stdin.write(b'1\r\n')
            # self.process.stdin.flush()

class Vector(object):
    def __init__(self, _point: list[float]):
        self.x = _point[0]
        self.y = _point[1]
        self.z = _point[2]

    def __str__(self):
        return "({}, {}, {})".format(self.x, self.y, self.z)

    def __iter__(self):
        return self

    def __next__(self):
        return self
    
    def radius(self):
        return (self.x**2 + self.y**2 + self.z**2)**0.5

class LightningTree(object):
    figure_tree = go.Figure()
    figure_plots = {}
    dimension = {'sum Q + q':'Кл', 'avg I':'А', 'avg Sigma':'См', 'sum Phi':'В', 'full Phi':'В'}
    eps0 = 8.85418781281313131313e-12
    k = 1 / (4 * math.pi * eps0)

    def __init__ (self, folder=None):
        if folder is None:
            return None
        self.folder = folder
        try:
            self.df_vertex = result = pd.read_csv(folder+'/vertex_table.txt', delim_whitespace=True)
            self.df_edge = pd.read_csv(folder+'/edge_table.txt', delim_whitespace=True)
            self.df_phi_info = pd.read_csv(folder+'/phi_info.txt', delim_whitespace=True)
        except pd.errors.EmptyDataError:
            print("Пустой датафрейм")
            self.default()
        except FileNotFoundError:
            print("Один из файлов не найден")
            self.default()
        if not self.df_edge.empty:
            self.df_edge['z'] = self.df_edge.apply(lambda row: (self.df_vertex[self.df_vertex['id'] == row['from']]['z'].item() + self.df_vertex[self.df_vertex['id'] == row['to']]['z'].item()) / 2, axis=1)
        else:
            self.df_edge['z'] = []

    def __str__(self):
        return "({}, {})".format(self.figure_tree, self.figure_plots)
    
    def default(self):
        self.df_vertex = pd.DataFrame([[1, 0, 0, 0, 0, 9000, 618974], [2,  0, 0, 0, 0, 9100, 1.99509e+07]], columns=["id", 'q', 'Q', 'x', 'y', 'z', "phi"])
        self.df_edge = pd.DataFrame([[1, 1, 2, 0.00060733, 1e-05]], columns=["id", "from", "to", "current", "sigma"])
        self.df_phi_info = pd.DataFrame([[7000, -5.80393e+07, -5.80393e+07], [7100, -6.5424e+07, -6.5424e+07]], columns=["z", "full_phi", "ext_phi"])
        self.plot_tree()
        self.plots()
        return self

    def distribution(self, df:pd.DataFrame, groupby:str, operarion:str, colomn:str) -> pd.DataFrame:
        """
        Метод расчета значений с группировкой данных по столбцу

        Parametrs
        ---------
        df: Данные для обработки
        groupby: Заголовок столбца по которому будет происходит группировка
        operation: Групповая операция по данным (sum, min, max, mean, ...)
        colomn: Заголовок столбца для которого будет осуществляться операция
        return: Сортированные по столбцу данные с распределением значений.
        """
        # return df.groupby([groupby]).agg({colomn:[operarion]})
        result = df.groupby([groupby]).agg({colomn:[operarion]})
        return pd.DataFrame({groupby : result.index.values[::-1], operarion+' '+colomn : result.values[::-1].ravel()})

    def fi_def(self, df:pd.DataFrame) -> pd.DataFrame:
        """
        Метод расчёта потенциала по столбцу
        
        Parametrs
        ---------
        df: Данные для обработки
        charge: заряд в вершине графа или в чехле ('q', 'Q')
        return: Сортированные по высоте данные с распределением потенциала
        """
        fi_list = []
        step = (int) ((df.z.max() - df.z.min()) / (df.nunique().z -1))
        for h in range(df.z.min(), df.z.max() + step, step):
            fi = 0
            for index, row in df.iterrows():
                if row.y == 0 and row.x == 0 and h == row.z:
                    fi += self.k * (row['q'] + row['Q']) * (1 / (Vector([0-row.x, 0 - row.y, h-row.z]).radius()+step/2)
                                        - 1 / (Vector([0-row.x, 0 - row.y, h+row.z]).radius()+step/2))
                else:
                    fi += self.k * (row['q'] + row['Q']) * (1 / (Vector([0-row.x, 0 - row.y, h-row.z]).radius())
                                        - 1 / (Vector([0-row.x, 0 - row.y, h+row.z]).radius()))
            fi_list.append([h, fi])
        
        return pd.DataFrame(fi_list, columns=['z', 'fi'])
        # fi = pd.DataFrame()          
        # for i, row in df.iterrows():
        #     fi[i] = df['z'].apply(lambda h: k * row[charge] * (1 / (Vector([0-row.x, 0 - row.y, h-row.z]).radius()+5) - 
        #                                                            1 / (Vector([0-row.x, 0 - row.y, h+row.z]).radius()+5))
        #                                         if row.y == 0 and row.x == 0 and h == row.z else
        #                                         k * row[charge] * (1 / (Vector([0-row.x, 0 - row.y, h-row.z]).radius()) - 
        #                                                            1 / (Vector([0-row.x, 0 - row.y, h+row.z]).radius()))
        #                             )
        # result = pd.DataFrame(df.z)
        # result['fi'] = fi.sum(axis=1)
        # return self.distribution(result,'z', 'mean', 'fi')

    def plot(self, _name:str, df:pd.DataFrame) -> go.Figure:
        """
        Создание 2D графика
        
        Parametrs
        ---------
        list_df: Лист данных для построения графика. (индекс - у, колонка - x)
        return: Графический объект
        """
        fig = go.Figure(go.Scatter(x=df[df.columns[1]], y=df[df.columns[0]], mode='lines+markers'),
                        layout=dict(uirevision=True, yaxis={'range':h_range}, xaxis={'title':_name + ', ' + self.dimension[_name]}))
        return fig

    def plots(self):
        """
        Создания словаря графических объектов из 2D графиков
        """
        self.figure_plots = {
            'sum Q + q': self.plot('sum Q + q', pd.DataFrame(self.distribution(self.df_vertex, 'z', 'sum', 'Q').merge(self.distribution(self.df_vertex, 'z', 'sum', 'q'), 'left').apply(lambda row: pd.Series([row[0], row[1] + row[2]], index=['z', 'sum Q + q']), axis=1))),
            'avg I': self.plot('avg I', self.distribution(self.df_edge, 'z', 'mean', 'current')),
            'avg Sigma': self.plot('avg Sigma', self.distribution(self.df_edge, 'z', 'mean', 'sigma')),
            'sum Phi': self.plot('sum Phi', self.fi_def(self.df_vertex)),
            'full Phi': self.plot('full Phi', self.distribution(self.df_phi_info, 'z', 'mean', 'full_phi'))

            # "sum_q" : self.plot([self.distribution(self.df_vertex, 'z', 'sum', 'q'), 
            #                     self.distribution(self.df_vertex, 'z', 'sum', 'Q')]),
            # "avg_q" : self.plot([self.distribution(self.df_vertex, 'z', 'mean', 'q'),
            #                          self.distribution(self.df_vertex, 'z', 'mean', 'Q')]),
            # 'Phi' : self.plot([self.distribution(self.df_phi_info, 'z', 'mean', 'full_phi'),
            #                    self.distribution(self.df_phi_info, 'z', 'mean', 'ext_phi')]),
            # 'current' : self.plot([self.distribution(self.df_edge, 'z', 'mean', 'current')]),
            # 'default' : self.plot([self.distribution(self.df_edge, 'z', 'mean', 'current'),
            #                        self.distribution(self.df_vertex, 'z', 'sum', 'Q'),
            #                        self.distribution(self.df_edge, 'z', 'mean', 'sigma'),
            #                        self.fi_def(self.df_vertex)
            #                        ])
            
        }

    def plot_tree(self):# -> go.Figure:
        """
        Создание графического объекта 3D графа дерева молнии
        """
        # Настройки для отображения графиков
        scale_nodes = [(0, "darkblue"), (0.15, "blue"), (0.49, "yellow"), (0.5, "gray"), (0.51, "yellow"), (0.85, "red"), (1, "darkred")] # цветовая шкала для зарядов
        scale_case = [(0, "darkblue"), (0.15, "blue"), (0.49, "yellow"), (0.5, "white"), (0.51, "yellow"), (0.85, "red"), (1, "darkred")] # цветовая шкала для чехлов
        setting = {'showbackground': False, 'showticklabels': True, 'showgrid': False, 'zeroline': True, 'range':[-2000, 2000]} # Параметры отображения системы координат
        setting_z = {'showbackground': True, 'showticklabels': True, 'showgrid': True, 'zeroline': True, 'range':h_range} # Параметры отображения системы координат для оси z

        # Создание настройки отображения графического объекта graph_object
        layout = go.Layout(showlegend=False, hovermode='closest',
                       scene={'xaxis': setting, 'yaxis': setting, 'zaxis': setting_z},
                       uirevision=True, scene_aspectmode='manual', scene_aspectratio=dict(x=1, y=1, z=2))
        
        # Набор DataFrame'а для создания графа
        array = []
        for index, row in self.df_edge.iterrows():
            array.append([row['from'], row.current])
            array.append([row['to'], row.current])
            array.append([None, None])
        df_edges = pd.DataFrame(array, columns=['id', 'current']).merge(self.df_vertex[['id', 'x', 'y', 'z']], on='id', how='left')
        # print(df_edges[df_edges.current.notna()].current.unique())
        
        # Построение рёбер
        edge_trace = go.Scatter3d(x=df_edges.x, y=df_edges.y, z=df_edges.z, 
                                  line=dict(width=2, color=df_edges.current, colorscale=["darkslateblue", "crimson"], cmin=-1, cmax=1),
                                  text=df_edges.current,
                                  hovertemplate='I=%{text}',
                                #   hoverinfo='none',
                                  mode='lines')
        
        # Построение зарядов    
        node_trace = go.Scatter3d(x=self.df_vertex.x, y=self.df_vertex.y, z=self.df_vertex.z,
                                  mode='markers',
                                  marker=dict(showscale=True, colorscale=scale_nodes, color=self.df_vertex.q, size=2.4),#, cmin=-0.001, cmax=0.001, ),
                                  line_width=.1)

        # Построение чехлов
        case_trace = go.Scatter3d(x=self.df_vertex.x, y=self.df_vertex.y, z=self.df_vertex.z,
                                  mode='markers',
                                  marker=dict(showscale=False, colorscale=scale_case, color=self.df_vertex.Q, size=12),#, cmin=-0.1, cmax=0.1),
                                  text = self.df_vertex.q,
                                  customdata= self.df_vertex.Q,
                                  hovertemplate='q= %{text} <br>Q= %{customdata}<extra></extra>',
                                  line_width=1,
                                  opacity=0.1)
        
        data = [edge_trace, node_trace, case_trace]
        self.figure_tree = go.Figure(data=data, layout=layout)
        # result = go.Figure(data=data, layout=layout)
        # return result


# # Создание Dash-приложения
app = Dash(__name__, server=server) #'SimpleExemple')
# Настройка и запуск Dash-приложения в зависимости от параметра
app.layout = html.Div([html.H1("Моделирование молнии", style={'textAlign': 'center', 'color': 'gold'}),
                    
                    dcc.Interval(id='interval-component', interval=1000, n_intervals=0, disabled=False),
                    
                    html.Div(["Параметр 1:", dcc.Input(id='param1', value='0.0', type='number'),
                                "Параметр 2:", dcc.Input(id='param2', value='0.0', type='number'),
                                "Параметр 3:", dcc.Input(id='param3', value='0.0', type='number'),
                                "Параметр 4: ",dcc.Input(id='param4', value='0.0', type='number'),
                                html.Button('Старт', id='start_button', n_clicks=0),
                                html.Button('Пауза', id='pause_button', n_clicks=0),
                                html.Button('Стоп', id='stop_button', n_clicks=0),
                                html.Button('Экспорт анимации', id='export_button', n_clicks=0)]),
                    html.Div(id="placeholder", style={"display":"none"}),
                    html.Div([html.Div([html.H4("Графики", style={'textAlign': 'center'}),
                                        html.Div([dcc.Graph(figure=LightningTree().default().figure_plots['sum Q + q'], id='plot_q', style={'height': '40vh', 'width': '33%'}),
                                                    dcc.Graph(figure=LightningTree().default().figure_plots['avg I'], id='plot_i', style={'height': '40vh', 'width': '32%'}),
                                                    dcc.Graph(figure=LightningTree().default().figure_plots['avg Sigma'], id='plot_sigma', style={'height': '40vh', 'width': '32%'})],
                                                    style={'display': 'flex', 'vertical-align':'top'}),
                                        html.Div([dcc.Graph(figure=LightningTree().default().figure_plots['sum Phi'], id='plot_phi', style={'height': '40vh', 'width': '49%'}),
                                                    dcc.Graph(figure=LightningTree().default().figure_plots['full Phi'], id='plot_full_phi', style={'height': '40vh', 'width': '49%'})],
                                                    style={'display': 'flex', 'vertical-align':'down'})],
                                        style={'width': '39%'}),

                                html.Div([html.H4("Граф дерева", style={'textAlign': 'center'}),
                                        dcc.Graph(figure=LightningTree().default().figure_tree, id='graph_tree', style={'height': '80vh'})],
                                        style={'width': '59%'})],
                                        style={'display':'flex'}),

                    html.Div(dcc.Slider(0, 1, step = 1, id='time_slider', disabled=False))], style={'width': '100%'})


class LightningTreePainter:

    def __init__(self, folder):
        self.lt_history = []
        self.folder = folder

    # @app.callback(Output('graph_tree', 'figure'))
    def update_graph(self):
        self.lt_history.append(LightningTree(self.folder))
        self.lt_history[-1].plot_tree()
        self.lt_history[-1].plots()
        return self.lt_history[-1].figure_tree, len(self.lt_history) - 1
    
    # @app.callback(Output('plot_q', 'figure'),
    #                 Output('plot_i', 'figure'),
    #                 Output('plot_sigma', 'figure'),
    #                 Output('plot_phi', 'figure'),
    #                 Output('plot_full_phi', 'figure'))
    def update_plots(self):
        return list(self.lt_history[-1].figure_plots.values())

    def subprocess_callback(self):
        self.update_graph()
        self.update_plots()
        update_layout()

class State(Enum):
        RUNNING = 1
        STOPPED = 2
        EXPORT = 3

class ActionController:

    def __init__(self, subprocess_controller, lightning_tree_painter):
        self.state = State.STOPPED
        self.subprocess_controller = subprocess_controller
        self.lightning_tree_painter = lightning_tree_painter


    def out_animations(self, lt_history:list[LightningTree], folder:str, format:str='jpg'):
        print('Запущен экпорт анимации')
        set_images = []
        w, h = 700, 500
        H = h*2
        for lt in lt_history:
            pic_tree = Image.open(io.BytesIO(pio.to_image(lt.figure_tree, 'png', width=w, height=H)))
            result = Image.new(pic_tree.mode, (w*2, H), color=0)

            for i, plot in enumerate(lt.figure_plots.values()):
                if i < 3:
                    pic_plots = Image.open(io.BytesIO(pio.to_image(plot, format, width=int(w/3))))
                    result.paste(pic_plots, box=(int(w/3)*i, 0))
                else:
                    pic_plots = Image.open(io.BytesIO(pio.to_image(plot, format, width=int(w/2))))
                    result.paste(pic_plots, box=(int(w/2)*(i-3), h))

            result.paste(pic_tree, box=(w, 0))

            set_images.append(result)
                
        set_images[0].save(
            folder + "/LightningTree.gif",
            save_all=True,
            append_images=set_images[1:],
            optimize=True,
            duration=500,
            loop=0)
        print('Экпорт завершён')
    
    
    def on_start_clicked(self):
        print("start_clicked")
        if self.state == State.EXPORT:
            print("Экспорт не завершен")
            return
        elif self.state == State.RUNNING:
            print("Процесс уже запущен")
        else:
            self.subprocess_controller.start(lightning_tree_painter.subprocess_callback)
            self.state = State.RUNNING

    def on_export_clicked(self):
        if self.state == State.EXPORT:
            print("Экспорт не завершен")
            return
        else:
            self.subprocess_controller.stop()
            self.state = State.EXPORT
            self.out_animations(lightning_tree_painter.lt_history, folder + "/Animations", 'jpg')
            self.state = State.STOPPED

    def on_stop_clicked(self):
        if self.state == State.EXPORT:
            print("Экспорт не завершен")
            return
        else:
            self.subprocess_controller.stop()
            self.state = State.STOPPED


subprocess_controller = SubprocessController()
lightning_tree_painter = LightningTreePainter( "LightningTree_data")
action_controller = ActionController(subprocess_controller, lightning_tree_painter)

@app.callback(Output(component_id='placeholder', component_property='children'),
              Input(component_id='start_button', component_property='n_clicks'))
def start(n_clicks):
    if n_clicks is not None:
        action_controller.on_start_clicked()
    return None

# @app.callback(Output(component_id='placeholder', component_property='children'),
#               Input(component_id='export_button', component_property='n_clicks'))
def export(n_clicks):
    if n_clicks is not None:
        action_controller.on_export_clicked()
    return None

# @app.callback(Output(component_id='placeholder', component_property='children'),
#               Input(component_id='stop_button', component_property='n_clicks'))
def stop(n_clicks):
    if n_clicks is not None:
        action_controller.on_stop_clicked()
    return None

def update_layout():
    app.layout = html.Div([html.H1("Моделирование молнии", style={'textAlign': 'center', 'color': 'gold'}),
                    
                    dcc.Interval(id='interval-component', interval=1000, n_intervals=0, disabled=False),
                    
                    html.Div(["Параметр 1:", dcc.Input(id='param1', value='0.0', type='number'),
                                "Параметр 2:", dcc.Input(id='param2', value='0.0', type='number'),
                                "Параметр 3:", dcc.Input(id='param3', value='0.0', type='number'),
                                "Параметр 4: ",dcc.Input(id='param4', value='0.0', type='number'),
                                html.Button('Старт', id='start_button', n_clicks=0),
                                html.Button('Пауза', id='pause_button', n_clicks=0),
                                html.Button('Стоп', id='stop_button', n_clicks=0),
                                html.Button('Экспорт анимации', id='export_button', n_clicks=0)]),
                    html.Div(id="placeholder", style={"display":"none"}),
                    html.Div([html.Div([html.H4("Графики", style={'textAlign': 'center'}),
                                        html.Div([dcc.Graph(figure=lightning_tree_painter.lt_history[-1].figure_plots['sum Q + q'], id='plot_q', style={'height': '40vh', 'width': '33%'}),
                                                    dcc.Graph(figure=lightning_tree_painter.lt_history[-1].figure_plots['avg I'], id='plot_i', style={'height': '40vh', 'width': '32%'}),
                                                    dcc.Graph(figure=lightning_tree_painter.lt_history[-1].figure_plots['avg Sigma'], id='plot_sigma', style={'height': '40vh', 'width': '32%'})],
                                                    style={'display': 'flex', 'vertical-align':'top'}),
                                        html.Div([dcc.Graph(figure=lightning_tree_painter.lt_history[-1].figure_plots['sum Phi'], id='plot_phi', style={'height': '40vh', 'width': '49%'}),
                                                    dcc.Graph(figure=lightning_tree_painter.lt_history[-1].figure_plots['full Phi'], id='plot_full_phi', style={'height': '40vh', 'width': '49%'})],
                                                    style={'display': 'flex', 'vertical-align':'down'})],
                                        style={'width': '39%'}),

                                html.Div([html.H4("Граф дерева", style={'textAlign': 'center'}),
                                        dcc.Graph(figure=lightning_tree_painter.lt_history[-1].figure_tree, id='graph_tree', style={'height': '80vh'})],
                                        style={'width': '59%'})],
                                        style={'display':'flex'}),

                    html.Div(dcc.Slider(0, 1, step = 1, id='time_slider', disabled=False))], style={'width': '100%'})

@server.route('/')
def run():
    """
    Метод для запуска Dash-приложения

    Parametrs
    ---------
    mode: Параметр запуска (inline - внутри jupyter; external - в браузере)
    interval: интервал обновления в секундах
    """
    app.run_server(debug = True)
    print("stopped")


if __name__ == '__main__':
    run()