#!/usr/bin/env python
# -*- coding: ansi_x3.4-1968 -*-
# generated by wxGlade 0.6.3 on Thu Aug 27 10:02:22 2009
import os
import wx
import chianti

# begin wxGlade: extracode
# end wxGlade



class ui_choice2Dialog(wx.Dialog):
    def __init__(self, *args, **kwds):
        # begin wxGlade: ui_choice2Dialog.__init__
        kwds["style"] = wx.DIALOG_NO_PARENT
        wx.Dialog.__init__(self, *args, **kwds)
        self.numListBox = wx.ListBox(self, -1, choices=[], style=wx.LB_MULTIPLE)
        self.button_ok = wx.Button(self, wx.ID_OK, "")
        self.denListBox = wx.ListBox(self, -1, choices=[], style=wx.LB_MULTIPLE)
        self.button_cancel = wx.Button(self, wx.ID_CANCEL, "")

        self.__set_properties()
        self.__do_layout()
        # end wxGlade

    def __set_properties(self):
        # begin wxGlade: ui_choice2Dialog.__set_properties
        self.SetTitle("ChiantiPy")
        _icon = wx.EmptyIcon()
        imagefile = os.path.join(chianti.__path__[0], "images/chianti2.png")
        _icon.CopyFromBitmap(wx.Bitmap(imagefile, wx.BITMAP_TYPE_ANY))
        self.SetIcon(_icon)
        self.SetSize((448, 556))
        # end wxGlade

    def __do_layout(self):
        # begin wxGlade: ui_choice2Dialog.__do_layout
        sizer_1 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_3 = wx.BoxSizer(wx.VERTICAL)
        sizer_2 = wx.BoxSizer(wx.VERTICAL)
        sizer_2.Add(self.numListBox, 1, wx.ALL|wx.EXPAND|wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 20)
        sizer_2.Add(self.button_ok, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 10)
        sizer_1.Add(sizer_2, 1, wx.EXPAND, 0)
        sizer_3.Add(self.denListBox, 1, wx.ALL|wx.EXPAND|wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 20)
        sizer_3.Add(self.button_cancel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 10)
        sizer_1.Add(sizer_3, 1, wx.EXPAND, 0)
        self.SetSizer(sizer_1)
        self.Layout()
        # end wxGlade

# end of class ui_choice2Dialog


#if __name__ == "__main__":
    #app = wx.PySimpleApp(0)
    #wx.InitAllImageHandlers()
    #ChiantiPy = ui_choice2Dialog(None, -1, "")
    #app.SetTopWindow(ChiantiPy)
    #ChiantiPy.Show()
    #app.MainLoop()
