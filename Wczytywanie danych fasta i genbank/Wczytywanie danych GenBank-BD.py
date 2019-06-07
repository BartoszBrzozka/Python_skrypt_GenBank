#!/usr/bin/env python
#  coding: cp1250 

import wx
from Bio import SeqIO
from Bio.Seq import Seq
import sys


def Zamknij(evt):
    dialog = wx.MessageDialog(Okno1, 'Czy na pewno?', 'Kończymy pracę', style = wx.OK | wx.CANCEL)
    x = dialog.ShowModal()
    dialog.Destroy()
    if x == wx.ID_OK:
        Okno1.Close()
def OtworzFasta(evt):
    dialog=wx.FileDialog(Okno1,message='Wybierz plik FASTA',defaultFile='',wildcard='*.FASTA',style=wx.FD_OPEN, pos=(10,10))
    if dialog.ShowModal() == wx.ID_OK:
        plik = dialog.GetPaths()
        print(plik)
        plik_tekst=open(plik[0],'r')
        Odczyt=plik_tekst.read()
        Odczyt=Odczyt.splitlines()
        plik_tekst.close()
        Sekwencje[ : ]=[]
        Nazwy[ : ]=[]
        nazwa=''
        for p in Odczyt:
            if p[0]=='>':
                if nazwa!='':
                    Sekwencje.append([gi,gb,nazwa,sekw])
                    Nazwy.append(nazwa)
                pp=p.split('|')
                gi=pp[2]
                if len(pp)>3:
                    gb=pp[1]
                    nazwa=pp[4]
                else:
                    gb=""
                    nazwa=pp[3]
                sekw=''
            else:
                sekw=sekw+p
        print (nazwa)
        Sekwencje.append([gi,gb,nazwa,sekw])
        Nazwy.append(nazwa)
        Lista.InsertItems(Nazwy,0)
        Lista.Show()
        dialog.Destroy()

def OtworzGenbank(evt):
    dialog=wx.FileDialog(Okno1,message='Wybierz plik GenBank',defaultFile='',wildcard='*.seq',style=wx.FD_OPEN, pos=(10,10))
    if dialog.ShowModal() == wx.ID_OK:
        plik = dialog.GetPath()
        for record in SeqIO.parse(plik, "genbank"):
            for y in record.features:
                if y.type == "CDS":
                    if "translation" in y.qualifiers:
                        Nazwy.append(record.id)
                        Desc.append(record.description)
                        Org.append(record.annotations['organism'])
                        if "product" in y.qualifiers:
                            Produkt.append(y.qualifiers['product'][0])
                        else:
                            Produkt.append('')
                        sekwencja = y.extract(record.seq)
                        if (len(sekwencja))%3==0:
                            protein_sequence = sekwencja.translate(table=11, cds=False, to_stop=True)
                            if(protein_sequence == y.qualifiers["translation"][0]):
                                d = 'Sekwencje zgodne'
                                Sekwencje.append(d)
                            else:
                                d = 'Sekwencje niezgodne'
                                Sekwencje.append(d)
                        elif (len(sekwencja)+1)%3==0:
                            sekwencja = sekwencja + Seq('N')
                            protein_sequence = sekwencja.translate(table=11, cds=False, to_stop=True)
                            if(protein_sequence == y.qualifiers["translation"][0]):
                                d = 'Sekwencje zgodne'
                                Sekwencje.append(d)
                            else:
                                d = 'Sekwencje niezgodne - brak nukleotydu'
                                Sekwencje.append(d)
                        elif (len(sekwencja)+2)%3==0:
                            sekwencja = sekwencja + Seq('NN')
                            protein_sequence = sekwencja.translate(table=11, cds=False, to_stop=True)
                            if(protein_sequence == y.qualifiers["translation"][0]):
                                d = 'Sekwencje zgodne'
                                Sekwencje.append(d)
                            else:
                                d = 'Sekwencje niezgodne - brak 2 nukleotydow'
                                Sekwencje.append(d)
                    else:
                        Nazwy.append(record.id)
                        Desc.append(record.description)
                        Org.append(record.annotations['organism'])
                        d = 'Pseudogen - brak sekwencji aminokwasowej'
                        Sekwencje.append(d)
                        if "product" in y.qualifiers:
                            Produkt.append(y.qualifiers['product'][0])
                        elif "note" in y.qualifiers:
                            Produkt.append(y.qualifiers['note'][0])
                        else:
                            Produkt.append('')
                elif y.type == "source":
                    if y.qualifiers["mol_type"][0] == "mRNA":
                        Nazwy.append(record.id)
                        Desc.append(record.description)
                        Org.append(record.annotations['organism'])
                        d = 'brak sekwencji aminokwasowej'
                        Sekwencje.append(d)
                        Produkt.append('')
                    else:
                        continue
            
        displaylist = list(zip(Nazwy,Desc,Org,Produkt,Sekwencje))
        list_ctrl.InsertColumn(0, 'Locus', width = 100)
        list_ctrl.InsertColumn(1, 'Opis', width = 400)
        list_ctrl.InsertColumn(2, 'Organizm', width = 200)
        list_ctrl.InsertColumn(3, 'Produkt', width = 150)
        list_ctrl.InsertColumn(4, 'Sekwencje', width = 150)
        index = 0
        for data in displaylist:
            list_ctrl.InsertItem(index , data[0])
            list_ctrl.SetItem(index, 1, data[1])
            list_ctrl.SetItem(index, 2, data[2])
            list_ctrl.SetItem(index, 3, data[3])
            list_ctrl.SetItem(index, 4, data[4])

        list_ctrl.Show()
        dialog.Destroy()
        

     
def IleNukleotydow(evt):
    Ktory=Lista.GetSelection()
    
    print(Ktory)
    print(Sekwencje[Ktory][3])
   
Prog=wx.App()
Produkt=[]
Sekwencje=[]
Nazwy=[]
displaylist = []
Desc = []
Org = []

Okno1=wx.Frame(None,title='Menu programu',size=(1300,600),pos=(200,50))

MenuListwa=wx.MenuBar()
ProgMenu=wx.Menu()
ProgMenuItem1=ProgMenu.Append(wx.ID_ANY,'Fasta','Czytaj Dane')
Okno1.Bind(wx.EVT_MENU, OtworzFasta, ProgMenuItem1)
ProgMenuItem2=ProgMenu.Append(wx.ID_ANY,'GenBank','Czytaj Dane')
Okno1.Bind(wx.EVT_MENU, OtworzGenbank, ProgMenuItem2)
MenuListwa.Append(ProgMenu,'Dane')

ProgMenu=wx.Menu()
ProgMenuItem1=ProgMenu.Append(wx.ID_ANY,'Liczba','Oblicz1')
Okno1.Bind(wx.EVT_MENU, IleNukleotydow, ProgMenuItem1)
MenuListwa.Append(ProgMenu,'Obliczenia')
ProgMenu=wx.Menu()
ProgMenuItem1=ProgMenu.Append(wx.ID_EXIT, 'Koniec', 'Koniec programu')
MenuListwa.Append(ProgMenu,'Wyjście')
Okno1.Bind(wx.EVT_MENU, Zamknij, ProgMenuItem1)

Okno1.SetMenuBar(MenuListwa)

panel = wx.Panel(parent = Okno1)
Lista=wx.ListBox(parent=panel, pos=(20,20), size=(400,400))
Lista1=wx.ListBox(parent=panel, pos=(230,20), size=(200,400))
list_ctrl=wx.ListCtrl(parent=panel, pos=(20,20), size=(1200,500), style = wx.LC_REPORT|wx.LC_HRULES|wx.LC_VRULES)
Lista.Hide()
Lista1.Hide()
list_ctrl.Hide()

Okno1.Show()


Prog.MainLoop()



