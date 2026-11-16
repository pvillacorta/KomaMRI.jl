module KomaMRSeqStudioExt

using Bonito
using Observables
import KomaMRI: KomaWindow, set_content!, _attach_sequence_studio!
import MRSeqStudio: Sequence, studio, StudioUI
using MRSeqStudio: StudioSequence

const STUDIOS = IdDict{KomaWindow,Observable{StudioSequence}}()
const PANELS = IdDict{KomaWindow,Any}()

function _attach_sequence_studio!(w::KomaWindow)
    sequence = Observable(StudioSequence())
    STUDIOS[w] = sequence
    PANELS[w] = StudioUI(sequence)
    w.handlers["view_sequence_studio"] = function (event; kwargs...)
        show_panel!(w)
        return nothing
    end
    push!(w.on_render, session -> Bonito.on_document_load(session, _menu_js()))
    _inject_menu!(w)
    push!(w.listeners, Observables.on(sequence) do studio_seq
        # Koma shows the sequence plot whenever w.seq changes. Put the same studio
        # node back if that assignment happened while the studio was open.
        restore = w.state[] == "sequence_studio"
        if !isempty(studio_seq.blocks)
            w.seq[] = Sequence(studio_seq)
        end
        restore && show_panel!(w)
    end)
    return nothing
end

function studio(w::KomaWindow)
    haskey(STUDIOS, w) || error("Open KomaUI with MRSeqStudio loaded before editing the studio.")
    return STUDIOS[w]
end

function show_panel!(w)
    set_content!(w, PANELS[w], "sequence_studio")
    return nothing
end

function _inject_menu!(w)
    session = w.session[]
    (isnothing(session) || !isopen(session)) && return nothing
    return Bonito.evaljs(session, _menu_js())
end

_menu_js() = js"""(() => {
    const mount = () => {
        if (document.getElementById('button_sequence_studio')) return;
        const defs = document.getElementById('button_sequence_definitions');
        if (!defs) {
            requestAnimationFrame(mount);
            return;
        }
        const item = document.createElement('li');
        const link = document.createElement('a');
        link.id = 'button_sequence_studio';
        link.className = 'dropdown-item';
        link.textContent = 'MR Seq Studio';
        link.href = '#';
        link.onclick = event => {
            event.preventDefault();
            KomaUI.notify('view_sequence_studio');
        };
        item.appendChild(link);
        defs.parentElement.insertAdjacentElement('afterend', item);
    };
    mount();
})()"""

end
