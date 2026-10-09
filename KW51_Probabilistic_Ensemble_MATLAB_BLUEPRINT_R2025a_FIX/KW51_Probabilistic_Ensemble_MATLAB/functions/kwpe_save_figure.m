function kwpe_save_figure(fig,base)
% Robust MATLAB R2025a-compatible figure export.
assert(isgraphics(fig,'figure'),'Invalid or deleted figure handle passed to kwpe_save_figure.');

% Force pending graphics operations to complete before export.
drawnow;

% Hide interactive axes toolbars so they are not embedded in exports.
axs=findall(fig,'Type','axes');
for k=1:numel(axs)
    try
        if isprop(axs(k),'Toolbar') && ~isempty(axs(k).Toolbar)
            axs(k).Toolbar.Visible='off';
        end
    catch
        % Toolbar availability varies across MATLAB graphics versions.
    end
end

drawnow;

pngFile=[base '.png'];
pdfFile=[base '.pdf'];

% exportgraphics is preferred over print for figures containing tiledlayout
% objects in current MATLAB releases (including R2025a).
try
    exportgraphics(fig,pngFile, ...
        'Resolution',220, ...
        'BackgroundColor','white');
catch MEpng
    warning('KW51:PNGExportFallback', ...
        'exportgraphics PNG failed (%s). Falling back to getframe/imwrite.',MEpng.message);
    frame=getframe(fig);
    imwrite(frame.cdata,pngFile);
end

try
    exportgraphics(fig,pdfFile, ...
        'ContentType','vector', ...
        'BackgroundColor','white');
catch MEvec
    warning('KW51:PDFVectorFallback', ...
        'Vector PDF export failed (%s). Retrying as image PDF.',MEvec.message);
    exportgraphics(fig,pdfFile, ...
        'ContentType','image', ...
        'Resolution',220, ...
        'BackgroundColor','white');
end
end
